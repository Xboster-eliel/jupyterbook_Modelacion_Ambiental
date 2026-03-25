import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.tsa.statespace.sarimax import SARIMAX


SARIMAX_CANDIDATES = [
    {"name": "ARIMAX(1,0,0)", "family": "ARIMAX", "order": (1, 0, 0), "seasonal_order": (0, 0, 0, 0)},
    {"name": "ARIMAX(1,0,1)", "family": "ARIMAX", "order": (1, 0, 1), "seasonal_order": (0, 0, 0, 0)},
    {"name": "ARIMAX(2,0,1)", "family": "ARIMAX", "order": (2, 0, 1), "seasonal_order": (0, 0, 0, 0)},
    {"name": "SARIMAX(1,0,0)x(1,0,0,12)", "family": "SARIMAX", "order": (1, 0, 0), "seasonal_order": (1, 0, 0, 12)},
    {"name": "SARIMAX(1,0,0)x(0,1,1,12)", "family": "SARIMAX", "order": (1, 0, 0), "seasonal_order": (0, 1, 1, 12)},
    {"name": "SARIMAX(1,0,1)x(0,1,1,12)", "family": "SARIMAX", "order": (1, 0, 1), "seasonal_order": (0, 1, 1, 12)},
    {"name": "SARIMAX(2,0,1)x(0,1,1,12)", "family": "SARIMAX", "order": (2, 0, 1), "seasonal_order": (0, 1, 1, 12)},
]


def read_precipitation(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=";")
    df.columns = [c.strip() for c in df.columns]
    for col in df.columns:
        if df[col].dtype == object:
            df[col] = df[col].astype(str).str.strip()
    df["date"] = pd.to_datetime(
        dict(
            year=pd.to_numeric(df["YYYY"], errors="coerce"),
            month=pd.to_numeric(df["MM"], errors="coerce"),
            day=pd.to_numeric(df["DD"], errors="coerce"),
        ),
        errors="coerce",
    )
    df["p_mm"] = pd.to_numeric(df["P [mm]"], errors="coerce")
    return df[["date", "p_mm"]].dropna().sort_values("date")


def read_flow(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df["date"] = pd.to_datetime(df["Fecha"], errors="coerce")
    df["q_m3s"] = pd.to_numeric(df["Valor"], errors="coerce")
    return df[["date", "q_m3s"]].dropna().sort_values("date")


def read_power_temperature(json_path: Path) -> pd.DataFrame:
    payload = json.loads(json_path.read_text(encoding="utf-8"))
    rows = []
    for key, value in payload["properties"]["parameter"]["T2M"].items():
        if key.endswith("13"):
            continue
        rows.append({"date": pd.to_datetime(key, format="%Y%m"), "t2m_c": float(value)})
    return pd.DataFrame(rows).sort_values("date")


def monthly_daylight_hours(lat_deg: float, month: int) -> float:
    lat_rad = math.radians(lat_deg)
    doy = pd.Timestamp(2001, month, 15).dayofyear
    delta = 0.409 * math.sin(2 * math.pi * doy / 365 - 1.39)
    omega = math.acos(max(-1.0, min(1.0, -math.tan(lat_rad) * math.tan(delta))))
    return 24 / math.pi * omega


def thornthwaite_pet(df_temp: pd.DataFrame, lat_deg: float) -> pd.DataFrame:
    df = df_temp.copy()
    heat_term = (df["t2m_c"].clip(lower=0) / 5) ** 1.514
    annual_i = heat_term.groupby(df.index.year).sum()
    annual_a = 6.75e-7 * annual_i**3 - 7.71e-5 * annual_i**2 + 1.792e-2 * annual_i + 0.49239
    daylight = {m: monthly_daylight_hours(lat_deg, m) for m in range(1, 13)}
    pet = []
    for date, row in df.iterrows():
        temp_c = max(row["t2m_c"], 0)
        heat = annual_i.loc[date.year]
        coeff = annual_a.loc[date.year]
        if temp_c <= 0 or heat == 0:
            pet.append(0.0)
        else:
            pet_mm = 16 * (daylight[date.month] / 12) * (date.days_in_month / 30) * ((10 * temp_c / heat) ** coeff)
            pet.append(float(pet_mm))
    df["pet_mm"] = pet
    return df


def prepare_monthly_dataframe(
    precip_daily: pd.DataFrame,
    flow_daily: pd.DataFrame,
    area_m2: float,
    temp_monthly: pd.DataFrame,
    lat_deg: float,
    date_start: str,
    date_end: str,
) -> pd.DataFrame:
    q = flow_daily.copy()
    q["q_mm_day"] = q["q_m3s"] * 86400.0 / area_m2 * 1000.0
    q_monthly = q.set_index("date").resample("MS").agg(
        q_mm_month=("q_mm_day", "sum"),
        q_mean_m3s=("q_m3s", "mean"),
        obs=("q_m3s", "count"),
    )
    q_monthly["days_in_month"] = q_monthly.index.days_in_month
    q_monthly["coverage"] = q_monthly["obs"] / q_monthly["days_in_month"]
    q_monthly.loc[q_monthly["coverage"] < 0.90, ["q_mm_month", "q_mean_m3s"]] = np.nan

    p_monthly = precip_daily.set_index("date").resample("MS").agg(p_mm=("p_mm", "sum"))
    t_monthly = thornthwaite_pet(temp_monthly.set_index("date").sort_index(), lat_deg)

    df = p_monthly.join(t_monthly[["t2m_c", "pet_mm"]], how="inner")
    df = df.join(q_monthly[["q_mm_month", "q_mean_m3s", "coverage"]], how="left")
    df["days_in_month"] = df.index.days_in_month
    df["q_mean_m3s_from_mm"] = df["q_mm_month"] * area_m2 / (1000.0 * 86400.0 * df["days_in_month"].astype(float))
    return df.loc[date_start:date_end].copy()


def train_test_split_time(df: pd.DataFrame, frac_train: float = 0.75):
    split_index = int(len(df) * frac_train)
    return df.iloc[:split_index].copy(), df.iloc[split_index:].copy(), split_index


def compute_metrics(obs, sim):
    obs = np.asarray(obs, dtype=float)
    sim = np.asarray(sim, dtype=float)
    mask = np.isfinite(obs) & np.isfinite(sim)
    obs = obs[mask]
    sim = sim[mask]
    rmse = float(np.sqrt(np.mean((sim - obs) ** 2)))
    mae = float(np.mean(np.abs(sim - obs)))
    bias = float(np.mean(sim - obs))
    nse = float(1 - np.sum((obs - sim) ** 2) / np.sum((obs - np.mean(obs)) ** 2))
    corr = float(np.corrcoef(obs, sim)[0, 1]) if len(obs) > 1 else np.nan
    return {"n": int(len(obs)), "rmse": rmse, "mae": mae, "bias": bias, "r2": nse, "nse": nse, "corr": corr}


def build_arx_design(df: pd.DataFrame, q_col: str = "q_mm_month") -> pd.DataFrame:
    out = df.copy()
    out["q_lag1"] = out[q_col].shift(1)
    out["q_lag2"] = out[q_col].shift(2)
    out["q_lag3"] = out[q_col].shift(3)
    out["sin_month"] = np.sin(2 * np.pi * out.index.month / 12)
    out["cos_month"] = np.cos(2 * np.pi * out.index.month / 12)
    return out


def fit_arx(train_df: pd.DataFrame):
    cols = ["q_lag1", "q_lag2", "q_lag3", "p_mm", "sin_month", "cos_month"]
    valid = train_df.dropna(subset=["q_mm_month"] + cols).copy()
    x = np.column_stack([np.ones(len(valid))] + [valid[c].to_numpy() for c in cols])
    y = valid["q_mm_month"].to_numpy()
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    return {"beta": beta, "columns": cols}


def predict_arx_in_sample(train_df: pd.DataFrame, model: dict) -> pd.DataFrame:
    cols = model["columns"]
    valid = train_df.dropna(subset=["q_mm_month"] + cols).copy()
    x = np.column_stack([np.ones(len(valid))] + [valid[c].to_numpy() for c in cols])
    valid["qhat_mm"] = np.maximum(x @ model["beta"], 0.0)
    return valid[["qhat_mm"]]


def predict_arx_recursive(df_full: pd.DataFrame, split_index: int, model: dict) -> pd.Series:
    q_work = df_full["q_mm_month"].copy()
    predictions = pd.Series(index=df_full.index, dtype=float)
    for i in range(split_index, len(df_full)):
        idx = df_full.index[i]
        vals = {
            "q_lag1": predictions.iloc[i - 1] if i - 1 >= split_index else q_work.iloc[i - 1],
            "q_lag2": predictions.iloc[i - 2] if i - 2 >= split_index else q_work.iloc[i - 2],
            "q_lag3": predictions.iloc[i - 3] if i - 3 >= split_index else q_work.iloc[i - 3],
            "p_mm": df_full.iloc[i]["p_mm"],
            "sin_month": math.sin(2 * math.pi * idx.month / 12),
            "cos_month": math.cos(2 * math.pi * idx.month / 12),
        }
        if not np.all(np.isfinite(list(vals.values()))):
            predictions.iloc[i] = np.nan
        else:
            x = np.r_[1.0, [vals[c] for c in model["columns"]]]
            predictions.iloc[i] = max(float(x @ model["beta"]), 0.0)
    return predictions


def run_water_balance(df_input: pd.DataFrame, k: float, s0: float = 0.0) -> pd.DataFrame:
    storage = float(s0)
    qsim = []
    stores = []
    for _, row in df_input.iterrows():
        available = max(storage + row["p_mm"] - row["pet_mm"], 0.0)
        q_t = min(k * available, available)
        storage = max(available - q_t, 0.0)
        qsim.append(float(q_t))
        stores.append(float(storage))
    out = df_input.copy()
    out["qsim_mm"] = qsim
    out["storage_mm"] = stores
    return out


def calibrate_water_balance(df_all: pd.DataFrame, split_index: int):
    best = None
    for k in np.linspace(0.05, 0.95, 181):
        sim = run_water_balance(df_all[["p_mm", "pet_mm"]], float(k))
        valid = df_all.iloc[:split_index][["q_mm_month"]].join(sim[["qsim_mm"]]).dropna()
        metrics = compute_metrics(valid["q_mm_month"], valid["qsim_mm"])
        if best is None or metrics["rmse"] < best["metrics"]["rmse"]:
            best = {"k": float(k), "metrics": metrics}
    return best


def fit_sarimax_candidates(train_all: pd.DataFrame, test_all: pd.DataFrame):
    rows = []
    best_by_family = {}
    for cand in SARIMAX_CANDIDATES:
        try:
            result = SARIMAX(
                train_all["q_mm_month"],
                exog=train_all[["p_mm"]],
                order=cand["order"],
                seasonal_order=cand["seasonal_order"],
                trend="c",
                enforce_stationarity=False,
                enforce_invertibility=False,
            ).fit(disp=False, maxiter=400)
            train_pred = result.fittedvalues
            test_pred = result.get_forecast(steps=len(test_all), exog=test_all[["p_mm"]]).predicted_mean
            row = {
                "name": cand["name"],
                "family": cand["family"],
                "aic": float(result.aic),
                "train_metrics": compute_metrics(train_all["q_mm_month"], train_pred),
                "test_metrics": compute_metrics(test_all["q_mm_month"], test_pred),
                "train_pred": train_pred,
                "test_pred": test_pred,
            }
            rows.append(row)
            current = best_by_family.get(cand["family"])
            if current is None or row["test_metrics"]["rmse"] < current["test_metrics"]["rmse"]:
                best_by_family[cand["family"]] = row
        except Exception as exc:
            rows.append({"name": cand["name"], "family": cand["family"], "aic": np.nan, "error": str(exc)})
    return rows, best_by_family


def summarize_case(
    precip_path: Path,
    flow_path: Path,
    power_json_path: Path,
    area_m2: float,
    lat_deg: float,
    date_start: str,
    date_end: str,
):
    precip_daily = read_precipitation(precip_path)
    flow_daily = read_flow(flow_path)
    temp_monthly = read_power_temperature(power_json_path)
    df = prepare_monthly_dataframe(precip_daily, flow_daily, area_m2, temp_monthly, lat_deg, date_start, date_end)
    train_all, test_all, split_index = train_test_split_time(df, frac_train=0.75)
    train_obs = train_all.dropna(subset=["q_mm_month"]).copy()
    test_obs = test_all.dropna(subset=["q_mm_month"]).copy()

    df_arx = build_arx_design(df)
    arx_model = fit_arx(df_arx.iloc[:split_index])
    df["q_arx_mm"] = predict_arx_recursive(df_arx, split_index, arx_model)
    train_fitted = predict_arx_in_sample(df_arx.iloc[:split_index], arx_model)
    df.loc[train_fitted.index, "q_arx_train_mm"] = train_fitted["qhat_mm"]
    df["q_arx_m3s"] = df["q_arx_mm"] * area_m2 / (1000.0 * 86400.0 * df["days_in_month"].astype(float))
    arx_train_metrics = compute_metrics(df.loc[train_fitted.index, "q_mm_month"], df.loc[train_fitted.index, "q_arx_train_mm"])
    arx_test_metrics = compute_metrics(test_obs["q_mm_month"], df.loc[test_obs.index, "q_arx_mm"])

    sarimax_rows, best_by_family = fit_sarimax_candidates(train_all, test_all)
    best_arimax = best_by_family["ARIMAX"]
    best_sarimax = best_by_family["SARIMAX"]
    df["q_arimax_mm"] = pd.concat([best_arimax["train_pred"], best_arimax["test_pred"]]).reindex(df.index)
    df["q_sarimax_mm"] = pd.concat([best_sarimax["train_pred"], best_sarimax["test_pred"]]).reindex(df.index)
    df["q_arimax_m3s"] = df["q_arimax_mm"] * area_m2 / (1000.0 * 86400.0 * df["days_in_month"].astype(float))
    df["q_sarimax_m3s"] = df["q_sarimax_mm"] * area_m2 / (1000.0 * 86400.0 * df["days_in_month"].astype(float))

    wb_best = calibrate_water_balance(df, split_index)
    wb_sim = run_water_balance(df[["p_mm", "pet_mm"]], wb_best["k"])
    df["q_phys_mm"] = wb_sim["qsim_mm"]
    df["storage_mm"] = wb_sim["storage_mm"]
    df["q_phys_m3s"] = df["q_phys_mm"] * area_m2 / (1000.0 * 86400.0 * df["days_in_month"].astype(float))
    phys_test_metrics = compute_metrics(test_obs["q_mm_month"], df.loc[test_obs.index, "q_phys_mm"])

    candidate_table = pd.DataFrame([
        {
            "modelo": row["name"],
            "familia": row["family"],
            "aic": row.get("aic", np.nan),
            "rmse_test": np.nan if row.get("test_metrics") is None else row["test_metrics"]["rmse"],
            "nse_test": np.nan if row.get("test_metrics") is None else row["test_metrics"]["nse"],
            "corr_test": np.nan if row.get("test_metrics") is None else row["test_metrics"]["corr"],
        }
        for row in sarimax_rows
    ]).sort_values(["familia", "rmse_test"])

    comparison = pd.DataFrame([
        {"modelo": "Estadistico ARX(3,1)", **arx_test_metrics},
        {"modelo": best_arimax["name"], **best_arimax["test_metrics"]},
        {"modelo": best_sarimax["name"], **best_sarimax["test_metrics"]},
        {"modelo": "Fisico balance + reservorio", **phys_test_metrics},
    ]).set_index("modelo").sort_values("rmse")

    return {
        "df": df,
        "train_all": train_all,
        "test_all": test_all,
        "train_obs": train_obs,
        "test_obs": test_obs,
        "split_index": split_index,
        "arx_model": arx_model,
        "arx_train_metrics": arx_train_metrics,
        "arx_test_metrics": arx_test_metrics,
        "best_arimax": best_arimax,
        "best_sarimax": best_sarimax,
        "wb_best": wb_best,
        "phys_test_metrics": phys_test_metrics,
        "candidate_table": candidate_table,
        "comparison": comparison,
    }


def render_text_table(df: pd.DataFrame) -> str:
    columns = list(df.columns)
    idx_name = df.index.name or "indice"
    idx_width = max(len(idx_name), max(len(str(idx)) for idx in df.index))
    widths = {}
    for col in columns:
        values = []
        for value in df[col]:
            if isinstance(value, (int, float, np.floating)) and pd.notna(value):
                values.append(f"{value:.3f}")
            else:
                values.append(str(value))
        widths[col] = max(len(col), max(len(v) for v in values))
    header = f"{idx_name.ljust(idx_width)} | " + " | ".join(col.ljust(widths[col]) for col in columns)
    sep = f"{'-' * idx_width}-|-" + "-|-".join("-" * widths[col] for col in columns)
    lines = [header, sep]
    for idx, row in df.iterrows():
        vals = []
        for col in columns:
            value = row[col]
            if isinstance(value, (int, float, np.floating)) and pd.notna(value):
                vals.append(f"{value:.3f}".ljust(widths[col]))
            else:
                vals.append(str(value).ljust(widths[col]))
        lines.append(f"{str(idx).ljust(idx_width)} | " + " | ".join(vals))
    return "\n".join(lines)


def fmt(value: float, digits: int = 2) -> str:
    return f"{value:.{digits}f}"
