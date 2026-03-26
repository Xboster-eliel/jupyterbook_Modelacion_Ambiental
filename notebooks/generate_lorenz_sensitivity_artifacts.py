from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import nbformat
from nbclient import NotebookClient
import numpy as np

from generate_lorenz_rk4_artifacts import BETA, PR, INITIAL, INITIAL_PERTURBED, integrate_lorenz

BASE_DIR = Path(__file__).resolve().parents[1]
FIG_DIR = BASE_DIR / 'content' / 'files_lorenz_sensibilidad'
SUMMARY_PATH = FIG_DIR / 'resumen_sensibilidad_lorenz.json'
NOTEBOOK_PATH = BASE_DIR / 'notebooks' / 'Guia_practica_Lorenz_Sensibilidad.ipynb'


def tail_mask(times: np.ndarray, start: float) -> np.ndarray:
    return times >= start


def sign_changes(series: np.ndarray) -> int:
    signs = np.sign(series)
    valid = signs != 0
    signs = signs[valid]
    if signs.size < 2:
        return 0
    return int(np.sum(signs[1:] != signs[:-1]))


def local_maxima(series: np.ndarray) -> np.ndarray:
    if series.size < 3:
        return np.array([], dtype=float)
    idx = (series[1:-1] > series[:-2]) & (series[1:-1] > series[2:])
    return series[1:-1][idx]


def occupancy_fraction(states: np.ndarray) -> float:
    return float(np.mean(states[:, 0] >= 0.0))


def summarize_window(times: np.ndarray, states: np.ndarray, start: float) -> dict[str, float]:
    mask = tail_mask(times, start)
    window = states[mask]
    return {
        'w_mean': float(window[:, 0].mean()),
        't1_mean': float(window[:, 1].mean()),
        't2_mean': float(window[:, 2].mean()),
        'w_std': float(window[:, 0].std()),
        't1_std': float(window[:, 1].std()),
        't2_std': float(window[:, 2].std()),
        't1_range': float(window[:, 1].max() - window[:, 1].min()),
        'occupancy_positive_w': occupancy_fraction(window),
        'sign_changes_w': sign_changes(window[:, 0]),
    }


def estimate_log_slope(times: np.ndarray, distances: np.ndarray, t_min: float, t_max: float) -> float:
    mask = (times >= t_min) & (times <= t_max) & (distances > 0)
    x = times[mask]
    y = np.log(distances[mask])
    if x.size < 2:
        return float('nan')
    slope, _ = np.polyfit(x, y, 1)
    return float(slope)


def resample_reference(ref_times: np.ndarray, ref_states: np.ndarray, target_times: np.ndarray) -> np.ndarray:
    interp = np.empty((target_times.size, 3), dtype=float)
    for i in range(3):
        interp[:, i] = np.interp(target_times, ref_times, ref_states[:, i])
    return interp


def make_ci_distance_plot(times: np.ndarray, states_a: np.ndarray, states_b: np.ndarray, path: Path) -> dict[str, float]:
    diff = np.linalg.norm(states_a - states_b, axis=1)
    eps = 1e-12
    fig, axes = plt.subplots(2, 1, figsize=(10.5, 7.5), sharex=True)
    axes[0].plot(times, diff, color='#111827', lw=1.2)
    for thresh, color in [(1.0, '#b91c1c'), (5.0, '#1d4ed8'), (10.0, '#0f766e')]:
        axes[0].axhline(thresh, color=color, linestyle='--', alpha=0.6, lw=0.9)
    axes[0].set_ylabel('||delta y||')
    axes[0].set_title('Separacion entre trayectorias para r = 30')
    axes[0].grid(alpha=0.25)

    axes[1].plot(times, np.log10(diff + eps), color='#7c3aed', lw=1.1)
    axes[1].set_xlabel('Tiempo')
    axes[1].set_ylabel('log10(||delta y||)')
    axes[1].grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)

    return {
        'time_gt_1': float(times[np.where(diff >= 1.0)[0][0]]),
        'time_gt_5': float(times[np.where(diff >= 5.0)[0][0]]),
        'time_gt_10': float(times[np.where(diff >= 10.0)[0][0]]),
        'log_slope_10_20': estimate_log_slope(times, diff, 10.0, 20.0),
        'log_slope_20_24': estimate_log_slope(times, diff, 20.0, 24.0),
        'final_distance': float(diff[-1]),
    }


def make_ci_component_plot(times: np.ndarray, states_a: np.ndarray, states_b: np.ndarray, path: Path) -> dict[str, float]:
    component_names = ['W', 'T1', 'T2']
    colors = ['#0f766e', '#b45309', '#1d4ed8']
    fig, axes = plt.subplots(3, 1, figsize=(10.5, 8.0), sharex=True)
    peaks = {}
    for i, (name, color) in enumerate(zip(component_names, colors)):
        delta = np.abs(states_a[:, i] - states_b[:, i])
        peaks[f'max_delta_{name.lower()}'] = float(delta.max())
        axes[i].plot(times, delta, color=color, lw=1.0)
        axes[i].set_ylabel(f'|delta {name}|')
        axes[i].grid(alpha=0.25)
        axes[i].set_title(f'Sensibilidad por componente: {name}', fontsize=10)
    axes[-1].set_xlabel('Tiempo')
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return peaks


def make_r_sweep_plots(r_values: np.ndarray, sweep: dict[float, dict[str, float]], metrics_path: Path, sign_path: Path, heatmap_path: Path) -> None:
    std_t1 = np.array([sweep[r]['t1_std'] for r in r_values])
    std_t2 = np.array([sweep[r]['t2_std'] for r in r_values])
    amp_t1 = np.array([sweep[r]['t1_range'] for r in r_values])
    occupancy = np.array([sweep[r]['occupancy_positive_w'] for r in r_values])
    sign_w = np.array([sweep[r]['sign_changes_w'] for r in r_values])
    sign_t1 = np.array([sweep[r]['sign_changes_t1'] for r in r_values])

    fig, axes = plt.subplots(2, 1, figsize=(10.5, 7.8), sharex=True)
    axes[0].plot(r_values, std_t1, color='#0f766e', lw=1.4, label='std(T1)')
    axes[0].plot(r_values, std_t2, color='#b45309', lw=1.4, label='std(T2)')
    axes[0].set_ylabel('Desviacion estandar')
    axes[0].legend(frameon=False)
    axes[0].grid(alpha=0.25)
    axes[1].plot(r_values, amp_t1, color='#1d4ed8', lw=1.4, label='rango(T1)')
    axes[1].plot(r_values, occupancy, color='#7c3aed', lw=1.2, label='fraccion W>=0')
    axes[1].set_xlabel('r')
    axes[1].set_ylabel('Metrica')
    axes[1].legend(frameon=False)
    axes[1].grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(metrics_path, dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    ax.plot(r_values, sign_w, color='#0f766e', lw=1.4, label='cambios de signo W')
    ax.plot(r_values, sign_t1, color='#b91c1c', lw=1.4, label='cambios de signo T1')
    ax.set_xlabel('r')
    ax.set_ylabel('Conteo')
    ax.set_title('Cambios de lobo y oscilacion segun r')
    ax.legend(frameon=False)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(sign_path, dpi=180)
    plt.close(fig)

    matrix = np.vstack([std_t1, std_t2, amp_t1, occupancy, sign_w]).astype(float)
    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    im = ax.imshow(matrix, aspect='auto', cmap='viridis', extent=[r_values.min(), r_values.max(), 4.5, -0.5])
    ax.set_yticks(range(5))
    ax.set_yticklabels(['std(T1)', 'std(T2)', 'rango(T1)', 'frac W>=0', 'signos W'])
    ax.set_xlabel('r')
    ax.set_title('Mapa de metricas para barrido fino en r')
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    fig.tight_layout()
    fig.savefig(heatmap_path, dpi=180)
    plt.close(fig)


def make_bifurcation_plot(r_values: np.ndarray, maxima: dict[float, np.ndarray], path: Path) -> None:
    fig, ax = plt.subplots(figsize=(10.8, 5.8))
    for r in r_values:
        vals = maxima[r]
        if vals.size:
            ax.plot(np.full(vals.shape, r), vals, '.', color='#111827', ms=1.6, alpha=0.7)
    ax.axvline(24.74, color='#b91c1c', linestyle='--', lw=1.0, alpha=0.7, label='umbral ~24.74')
    ax.set_xlabel('r')
    ax.set_ylabel('Maximos locales de T1')
    ax.set_title('Barrido tipo bifurcacion en r')
    ax.legend(frameon=False)
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def make_dt_plots(dt_results: dict[float, dict[str, np.ndarray | float]], short_path: Path, error_path: Path, stats_path: Path, phase_path: Path) -> None:
    dts = sorted(dt_results.keys(), reverse=True)
    colors = ['#b91c1c', '#0f766e', '#1d4ed8', '#7c3aed', '#111827']

    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    for dt, color in zip(dts, colors):
        times = dt_results[dt]['times']
        states = dt_results[dt]['states']
        mask = times <= 10.0
        ax.plot(times[mask], states[mask, 1], lw=1.0, color=color, label=f'dt={dt:g}')
    ax.set_xlabel('Tiempo')
    ax.set_ylabel('T1')
    ax.set_title('Sensibilidad de T1(t) al paso temporal')
    ax.legend(frameon=False, ncol=3)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(short_path, dpi=180)
    plt.close(fig)

    dt_vals = np.array(sorted(dt_results.keys()))
    l2 = np.array([dt_results[dt]['l2_error_t1'] for dt in dt_vals])
    linf = np.array([dt_results[dt]['linf_error_t1'] for dt in dt_vals])
    fig, ax = plt.subplots(figsize=(8.6, 4.8))
    ax.loglog(dt_vals, l2, 'o-', color='#0f766e', label='L2(T1)')
    ax.loglog(dt_vals, linf, 's-', color='#b91c1c', label='Linf(T1)')
    ax.invert_xaxis()
    ax.set_xlabel('dt')
    ax.set_ylabel('Error respecto a referencia')
    ax.set_title('Convergencia numerica frente a dt')
    ax.legend(frameon=False)
    ax.grid(alpha=0.25, which='both')
    fig.tight_layout()
    fig.savefig(error_path, dpi=180)
    plt.close(fig)

    labels = [f'{dt:g}' for dt in dt_vals]
    mean_t1 = [dt_results[dt]['mean_t1_last20'] for dt in dt_vals]
    std_t1 = [dt_results[dt]['std_t1_last20'] for dt in dt_vals]
    mean_t2 = [dt_results[dt]['mean_t2_last20'] for dt in dt_vals]
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.5))
    axes[0].plot(labels, mean_t1, 'o-', color='#1d4ed8', label='media T1')
    axes[0].plot(labels, std_t1, 's-', color='#0f766e', label='std T1')
    axes[0].legend(frameon=False)
    axes[0].grid(alpha=0.25)
    axes[0].set_title('Estadisticos de T1 en [20, 50]')
    axes[0].set_xlabel('dt')
    axes[1].plot(labels, mean_t2, 'o-', color='#b45309')
    axes[1].set_title('Media de T2 en [20, 50]')
    axes[1].set_xlabel('dt')
    axes[1].grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(stats_path, dpi=180)
    plt.close(fig)

    phase_dts = [0.02, 0.005, 0.0025]
    fig, axes = plt.subplots(1, 3, figsize=(14.2, 4.2), sharex=True, sharey=True)
    for ax, dt, color in zip(axes, phase_dts, ['#b91c1c', '#1d4ed8', '#111827']):
        states = dt_results[dt]['states']
        mask = dt_results[dt]['times'] >= 20.0
        ax.plot(states[mask, 0], states[mask, 1], color=color, lw=0.35)
        ax.set_title(f'dt = {dt:g}')
        ax.set_xlabel('W')
        ax.grid(alpha=0.25)
    axes[0].set_ylabel('T1')
    fig.suptitle('Comparacion cualitativa del atractor segun dt', y=1.02)
    fig.tight_layout()
    fig.savefig(phase_path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def make_transient_plots(transient_summary: dict[str, dict[str, dict[str, float]]], metrics_path: Path, series_path: Path, boxplot_path: Path) -> None:
    cases = ['24.0', '25.0', '30.0']
    windows = ['0', '10', '20']
    labels = ['[0,50]', '[10,50]', '[20,50]']

    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.6), sharey=True)
    for ax, case in zip(axes, cases):
        t1_std = [transient_summary[case][w]['t1_std'] for w in windows]
        t1_range = [transient_summary[case][w]['t1_range'] for w in windows]
        ax.plot(labels, t1_std, 'o-', color='#0f766e', label='std(T1)')
        ax.plot(labels, t1_range, 's-', color='#1d4ed8', label='rango(T1)')
        ax.set_title(f'r = {case[:-2]}')
        ax.grid(alpha=0.25)
        if ax is axes[0]:
            ax.set_ylabel('Metrica')
            ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(metrics_path, dpi=180)
    plt.close(fig)

    case = 24.0
    times, states = integrate_lorenz(case, INITIAL)
    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    ax.plot(times, states[:, 1], color='#1d4ed8', lw=0.8)
    ax.axvspan(0, 10, color='#f59e0b', alpha=0.15, label='transitorio temprano')
    ax.axvspan(10, 20, color='#10b981', alpha=0.12, label='zona intermedia')
    ax.axvspan(20, 50, color='#6366f1', alpha=0.10, label='ventana final')
    ax.set_xlabel('Tiempo')
    ax.set_ylabel('T1')
    ax.set_title('Importancia del transitorio en r = 24')
    ax.legend(frameon=False, ncol=3, fontsize=8)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(series_path, dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9.2, 4.6))
    data = [[transient_summary[c][w]['t1_std'] for c in cases] for w in windows]
    ax.boxplot(data, labels=labels)
    ax.set_ylabel('std(T1) por caso')
    ax.set_title('Comparacion de dispersion al descartar transitorio')
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(boxplot_path, dpi=180)
    plt.close(fig)


def make_observable_plots(observable_summary: dict[str, dict[str, float]], summary_path: Path, hist_path: Path, occupancy_path: Path) -> None:
    cases = ['24.0', '25.0', '30.0']
    labels = ['24', '25', '30']
    mean_t1 = [observable_summary[c]['t1_mean'] for c in cases]
    std_t1 = [observable_summary[c]['t1_std'] for c in cases]
    range_t1 = [observable_summary[c]['t1_range'] for c in cases]
    occupancy = [observable_summary[c]['occupancy_positive_w'] for c in cases]

    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.8))
    axes[0].plot(labels, mean_t1, 'o-', color='#1d4ed8', label='media T1')
    axes[0].plot(labels, std_t1, 's-', color='#0f766e', label='std T1')
    axes[0].plot(labels, range_t1, '^-', color='#b45309', label='rango T1')
    axes[0].legend(frameon=False, fontsize=8)
    axes[0].set_title('Observables resumidos por caso')
    axes[0].grid(alpha=0.25)
    axes[1].bar(labels, occupancy, color=['#0f766e', '#1d4ed8', '#7c3aed'])
    axes[1].set_ylim(0, 1)
    axes[1].set_title('Fraccion de tiempo con W >= 0')
    axes[1].grid(alpha=0.25, axis='y')
    fig.tight_layout()
    fig.savefig(summary_path, dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.4), sharey=True)
    colors = ['#0f766e', '#1d4ed8', '#7c3aed']
    for ax, case, color in zip(axes, [24.0, 25.0, 30.0], colors):
        times, states = integrate_lorenz(case, INITIAL)
        mask = times >= 20.0
        ax.hist(states[mask, 1], bins=35, density=True, color=color, alpha=0.75)
        ax.set_title(f'r = {int(case)}')
        ax.set_xlabel('T1')
        ax.grid(alpha=0.2)
    axes[0].set_ylabel('Densidad')
    fig.tight_layout()
    fig.savefig(hist_path, dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2), sharex=True, sharey=True, constrained_layout=True)
    for ax, case in zip(axes, [24.0, 25.0, 30.0]):
        times, states = integrate_lorenz(case, INITIAL)
        mask = times >= 20.0
        h = ax.hist2d(states[mask, 0], states[mask, 1], bins=50, cmap='magma')
        ax.set_title(f'r = {int(case)}')
        ax.set_xlabel('W')
    axes[0].set_ylabel('T1')
    fig.colorbar(h[3], ax=axes, fraction=0.025, pad=0.02)
    fig.savefig(occupancy_path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def image_cell(image_name: str) -> nbformat.NotebookNode:
    return nbformat.v4.new_code_cell(
        "image_path = BASE_DIR / 'content' / 'files_lorenz_sensibilidad' / '" + image_name + "'\n"
        'img = plt.imread(image_path)\n'
        'plt.figure(figsize=(11, 5.5))\n'
        'plt.imshow(img)\n'
        "plt.axis('off')\n"
        'plt.show()'
    )


def markdown_image_cell(image_path: str, alt_text: str) -> nbformat.NotebookNode:
    return nbformat.v4.new_markdown_cell(f'![{alt_text}]({image_path})')


def build_notebook() -> nbformat.NotebookNode:
    nb = nbformat.v4.new_notebook()
    nb.metadata['kernelspec'] = {'display_name': 'Python 3', 'language': 'python', 'name': 'python3'}
    nb.metadata['language_info'] = {'name': 'python', 'version': '3'}
    cells = []
    cells.append(nbformat.v4.new_markdown_cell(
        '# Guia practica: analisis de sensibilidad del sistema de Lorenz\n\n'
        'Este notebook reune el analisis ampliado de sensibilidad para el sistema de Lorenz con `RK4`, '
        'incluyendo sensibilidad a condiciones iniciales, barrido fino en `r`, sensibilidad a `dt`, '
        'transitorio y observables resumidos.'
    ))
    cells.append(nbformat.v4.new_code_cell(
        'import json\nfrom pathlib import Path\nimport matplotlib.pyplot as plt\n\n'
        "BASE_DIR = Path.cwd().resolve().parents[0]\n"
        "SUMMARY_PATH = BASE_DIR / 'content' / 'files_lorenz_sensibilidad' / 'resumen_sensibilidad_lorenz.json'\n"
        "with SUMMARY_PATH.open('r', encoding='utf-8') as f:\n    summary = json.load(f)\n\nsummary.keys()"
    ))
    sections = [
        ('## Sensibilidad a condiciones iniciales', 'lorenz_ci_distancia.png', 'La separacion entre trayectorias muestra como un cambio extremadamente pequeno en `T2(0)` produce una divergencia apreciable despues de un horizonte finito.'),
        ('## Sensibilidad por componente', 'lorenz_ci_componentes.png', 'El crecimiento no es uniforme en las tres variables. Esta desagregacion ayuda a ver en que componente se vuelve mas visible la perdida de sincronizacion.'),
        ('## Barrido fino en r', 'lorenz_r_barrido_metricas.png', 'Las metricas de dispersion y amplitud confirman que la vecindad `r = 24-25` es especialmente sensible.'),
        ('## Cambios de signo y cambios de lobo', 'lorenz_r_cambios_signo.png', 'El conteo de cambios de signo resume la frecuencia con la que la trayectoria cambia de lobo o de orientacion dinamica.'),
        ('## Mapa de metricas en r', 'lorenz_r_heatmap_metricas.png', 'La vista matricial permite identificar regiones de bajo y alto dinamismo sin depender de una sola metrica.'),
        ('## Barrido tipo bifurcacion', 'lorenz_bifurcacion_t1_vs_r.png', 'El diagrama de bifurcacion ayuda a distinguir zonas de equilibrio, transicion y comportamiento caotico.'),
        ('## Sensibilidad a dt en ventana corta', 'lorenz_dt_t1_corto.png', 'La convergencia puntual es razonable en tiempos cortos cuando se refina `dt`.'),
        ('## Error numerico frente a referencia', 'lorenz_dt_error_vs_dt.png', 'La disminucion del error al refinar `dt` da una medida concreta de robustez numerica.'),
        ('## Estadisticos segun dt', 'lorenz_dt_estadisticos.png', 'En sistemas caoticos, la consistencia estadistica es mas informativa que la coincidencia puntual en tiempos largos.'),
        ('## Comparacion del atractor segun dt', 'lorenz_dt_fases_comparacion.png', 'La forma global del atractor permanece reconocible cuando `dt` es suficientemente pequeno.'),
        ('## Sensibilidad al transitorio', 'lorenz_transitorio_metricas.png', 'Las metricas cambian de manera importante cuando se descarta el arranque, especialmente cerca del umbral.'),
        ('## Ventanas temporales y transitorio', 'lorenz_transitorio_series.png', 'Esta figura ayuda a ver por que una solucion puede parecer casi estacionaria en una ventana y no serlo globalmente.'),
        ('## Dispersion segun ventana', 'lorenz_transitorio_boxplot.png', 'Comparar ventanas evita confundir un transitorio largo con un equilibrio real.'),
        ('## Observables resumidos', 'lorenz_observables_resumidos.png', 'Las medias, dispersiones y fracciones de ocupacion permiten separar trayectoria puntual y estructura global.'),
        ('## Histogramas de T1', 'lorenz_histogramas_t1.png', 'La distribucion de `T1` en la ventana final muestra como cambia la ocupacion del atractor entre casos.'),
        ('## Densidad en el plano W-T1', 'lorenz_ocupacion_w_t1.png', 'Los mapas de densidad ilustran la ocupacion estadistica del atractor mas alla de una sola trayectoria.')
    ]
    for title, image, text in sections:
        cells.append(nbformat.v4.new_markdown_cell(title))
        cells.append(image_cell(image))
        cells.append(nbformat.v4.new_markdown_cell(text))
    cells.append(nbformat.v4.new_markdown_cell('## Tabla resumida'))
    cells.append(nbformat.v4.new_code_cell("summary['initial_condition_sensitivity']"))
    cells.append(nbformat.v4.new_code_cell("summary['dt_sensitivity']"))
    nb.cells = cells
    return nb


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    summary: dict[str, object] = {
        'parameters': {
            'Pr': PR,
            'b': BETA,
            'initial_condition': INITIAL.tolist(),
            'perturbed_initial_condition': INITIAL_PERTURBED.tolist(),
        }
    }

    # 1. Sensibilidad a condiciones iniciales
    times_30, states_30 = integrate_lorenz(30.0, INITIAL)
    _, states_30_b = integrate_lorenz(30.0, INITIAL_PERTURBED)
    ci_metrics = make_ci_distance_plot(times_30, states_30, states_30_b, FIG_DIR / 'lorenz_ci_distancia.png')
    ci_metrics.update(make_ci_component_plot(times_30, states_30, states_30_b, FIG_DIR / 'lorenz_ci_componentes.png'))
    summary['initial_condition_sensitivity'] = ci_metrics

    # 2. Sensibilidad paramétrica local en r
    r_values = np.round(np.arange(22.0, 28.01, 0.25), 2)
    sweep: dict[float, dict[str, float]] = {}
    for r in r_values:
        times, states = integrate_lorenz(float(r), INITIAL)
        mask = times >= 20.0
        window = states[mask]
        sweep[float(r)] = {
            't1_std': float(window[:, 1].std()),
            't2_std': float(window[:, 2].std()),
            't1_range': float(window[:, 1].max() - window[:, 1].min()),
            'occupancy_positive_w': occupancy_fraction(window),
            'sign_changes_w': sign_changes(window[:, 0]),
            'sign_changes_t1': sign_changes(window[:, 1]),
        }
    make_r_sweep_plots(
        r_values,
        sweep,
        FIG_DIR / 'lorenz_r_barrido_metricas.png',
        FIG_DIR / 'lorenz_r_cambios_signo.png',
        FIG_DIR / 'lorenz_r_heatmap_metricas.png',
    )
    summary['r_local_sensitivity'] = {str(k): v for k, v in sweep.items()}

    # 3. Barrido tipo bifurcacion
    r_bif = np.round(np.linspace(0.0, 40.0, 121), 3)
    bif_maxima: dict[float, np.ndarray] = {}
    for r in r_bif:
        times, states = integrate_lorenz(float(r), INITIAL, dt=0.01, t_stop=40.0)
        mask = times >= 20.0
        maxima = local_maxima(states[mask, 1])
        bif_maxima[float(r)] = maxima[:18]
    make_bifurcation_plot(r_bif, bif_maxima, FIG_DIR / 'lorenz_bifurcacion_t1_vs_r.png')
    summary['bifurcation'] = {
        'r_min': float(r_bif.min()),
        'r_max': float(r_bif.max()),
        'samples': int(r_bif.size),
    }

    # 4. Sensibilidad numerica a dt
    dt_values = [0.02, 0.01, 0.005, 0.0025]
    ref_times, ref_states = integrate_lorenz(30.0, INITIAL, dt=0.00125, t_stop=50.0)
    dt_results: dict[float, dict[str, np.ndarray | float]] = {}
    for dt in dt_values:
        times, states = integrate_lorenz(30.0, INITIAL, dt=dt, t_stop=50.0)
        ref_interp = resample_reference(ref_times, ref_states, times)
        mask_short = times <= 10.0
        diff_short = states[mask_short, 1] - ref_interp[mask_short, 1]
        mask_tail = times >= 20.0
        dt_results[dt] = {
            'times': times,
            'states': states,
            'l2_error_t1': float(np.sqrt(np.mean(diff_short ** 2))),
            'linf_error_t1': float(np.max(np.abs(diff_short))),
            'mean_t1_last20': float(states[mask_tail, 1].mean()),
            'std_t1_last20': float(states[mask_tail, 1].std()),
            'mean_t2_last20': float(states[mask_tail, 2].mean()),
            'std_t2_last20': float(states[mask_tail, 2].std()),
        }
    make_dt_plots(
        dt_results,
        FIG_DIR / 'lorenz_dt_t1_corto.png',
        FIG_DIR / 'lorenz_dt_error_vs_dt.png',
        FIG_DIR / 'lorenz_dt_estadisticos.png',
        FIG_DIR / 'lorenz_dt_fases_comparacion.png',
    )
    summary['dt_sensitivity'] = {
        str(dt): {
            'l2_error_t1': float(dt_results[dt]['l2_error_t1']),
            'linf_error_t1': float(dt_results[dt]['linf_error_t1']),
            'mean_t1_last20': float(dt_results[dt]['mean_t1_last20']),
            'std_t1_last20': float(dt_results[dt]['std_t1_last20']),
        }
        for dt in dt_values
    }

    # 5. Sensibilidad al transitorio
    transient_cases = [24.0, 25.0, 30.0]
    transient_summary: dict[str, dict[str, dict[str, float]]] = {}
    for case in transient_cases:
        times, states = integrate_lorenz(case, INITIAL)
        transient_summary[str(case)] = {
            '0': summarize_window(times, states, 0.0),
            '10': summarize_window(times, states, 10.0),
            '20': summarize_window(times, states, 20.0),
        }
    make_transient_plots(
        transient_summary,
        FIG_DIR / 'lorenz_transitorio_metricas.png',
        FIG_DIR / 'lorenz_transitorio_series.png',
        FIG_DIR / 'lorenz_transitorio_boxplot.png',
    )
    summary['transient_sensitivity'] = transient_summary

    # 6. Observables resumidos
    observable_summary: dict[str, dict[str, float]] = {}
    for case in transient_cases:
        times, states = integrate_lorenz(case, INITIAL)
        mask = times >= 20.0
        observable_summary[str(case)] = {
            't1_mean': float(states[mask, 1].mean()),
            't1_std': float(states[mask, 1].std()),
            't1_range': float(states[mask, 1].max() - states[mask, 1].min()),
            'occupancy_positive_w': occupancy_fraction(states[mask]),
        }
    make_observable_plots(
        observable_summary,
        FIG_DIR / 'lorenz_observables_resumidos.png',
        FIG_DIR / 'lorenz_histogramas_t1.png',
        FIG_DIR / 'lorenz_ocupacion_w_t1.png',
    )
    summary['observable_summary'] = observable_summary

    with SUMMARY_PATH.open('w', encoding='utf-8') as f:
        json.dump(summary, f, indent=2)

    nb = build_notebook()
    with NOTEBOOK_PATH.open('w', encoding='utf-8') as f:
        nbformat.write(nb, f)
    client = NotebookClient(nb, timeout=240, kernel_name='python3')
    client.execute(cwd=str(NOTEBOOK_PATH.parent))
    with NOTEBOOK_PATH.open('w', encoding='utf-8') as f:
        nbformat.write(nb, f)

    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
