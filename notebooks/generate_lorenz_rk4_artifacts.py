from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
import nbformat
from nbclient import NotebookClient
import numpy as np

SIGMA = 10.0
BETA = 8.0 / 3.0
DT = 0.005
T_STOP = 50.0
INITIAL = np.array([1.0, 1.0, 1.0], dtype=float)
INITIAL_PERTURBED = np.array([1.001, 1.0, 1.0], dtype=float)
BASE_DIR = Path(__file__).resolve().parents[1]
CONTENT_DIR = BASE_DIR / 'content'
FIG_DIR = CONTENT_DIR / 'files_lorenz_rk4'
NOTEBOOK_PATH = BASE_DIR / 'notebooks' / 'Guia_practica_Lorenz_RK4.ipynb'
SUMMARY_PATH = FIG_DIR / 'resumen_lorenz_rk4.json'


def lorenz_rhs(state: np.ndarray, r: float, sigma: float = SIGMA, beta: float = BETA) -> np.ndarray:
    x, y, z = state
    return np.array(
        [
            sigma * (y - x),
            x * (r - z) - y,
            x * y - beta * z,
        ],
        dtype=float,
    )


def rk4_step(state: np.ndarray, dt: float, r: float) -> np.ndarray:
    k1 = lorenz_rhs(state, r)
    k2 = lorenz_rhs(state + 0.5 * dt * k1, r)
    k3 = lorenz_rhs(state + 0.5 * dt * k2, r)
    k4 = lorenz_rhs(state + dt * k3, r)
    return state + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def integrate_lorenz(r: float, y0: np.ndarray, dt: float = DT, t_stop: float = T_STOP) -> tuple[np.ndarray, np.ndarray]:
    times = np.arange(0.0, t_stop + dt, dt)
    states = np.empty((times.size, 3), dtype=float)
    states[0] = y0
    for idx in range(times.size - 1):
        states[idx + 1] = rk4_step(states[idx], dt, r)
    return times, states


def fixed_points(r: float, beta: float = BETA) -> list[np.ndarray]:
    if r <= 1.0:
        return [np.array([0.0, 0.0, 0.0])]
    coord = np.sqrt(beta * (r - 1.0))
    return [
        np.array([0.0, 0.0, 0.0]),
        np.array([coord, coord, r - 1.0]),
        np.array([-coord, -coord, r - 1.0]),
    ]


def min_distance_to_nonzero_fixed_point(states: np.ndarray, r: float) -> float | None:
    points = fixed_points(r)
    if len(points) < 3:
        return None
    distances = [np.linalg.norm(states - point, axis=1).min() for point in points[1:]]
    return float(min(distances))


def divergence_time(states_a: np.ndarray, states_b: np.ndarray, times: np.ndarray, threshold: float = 5.0) -> float | None:
    distances = np.linalg.norm(states_a - states_b, axis=1)
    hits = np.where(distances >= threshold)[0]
    if hits.size == 0:
        return None
    return float(times[hits[0]])


def settling_window_std(series: np.ndarray, window: int = 2000) -> float:
    return float(np.std(series[-window:]))


def make_time_plot(times: np.ndarray, states: np.ndarray, r: float, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.plot(times, states[:, 0], color='#0d3b66', lw=1.0)
    ax.set_title(f'Evolucion temporal de x(t) para r = {r}')
    ax.set_xlabel('Tiempo')
    ax.set_ylabel('x(t)')
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def make_phase_plot(states: np.ndarray, r: float, path: Path) -> None:
    fig = plt.figure(figsize=(7, 5.5))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(states[:, 0], states[:, 1], states[:, 2], color='#7c3aed', lw=0.55)
    ax.set_title(f'Trayectoria en espacio de fases para r = {r}')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def make_comparison_time_plot(cases: dict[float, tuple[np.ndarray, np.ndarray]], path: Path) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    for ax, r, color in zip(axes, [10.0, 24.0], ['#0f766e', '#b45309']):
        times, states = cases[r]
        ax.plot(times, states[:, 0], color=color, lw=1.0)
        ax.set_title(f'x(t) para r = {int(r)}')
        ax.set_ylabel('x(t)')
        ax.grid(alpha=0.25)
    axes[-1].set_xlabel('Tiempo')
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def make_comparison_phase_plot(cases: dict[float, tuple[np.ndarray, np.ndarray]], path: Path) -> None:
    fig = plt.figure(figsize=(12, 5.5))
    for i, (r, color) in enumerate([(10.0, '#0f766e'), (24.0, '#b45309')], start=1):
        ax = fig.add_subplot(1, 2, i, projection='3d')
        states = cases[r][1]
        ax.plot(states[:, 0], states[:, 1], states[:, 2], color=color, lw=0.5)
        ax.set_title(f'Fases para r = {int(r)}')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def make_r25_combo_plot(times: np.ndarray, states: np.ndarray, path: Path) -> None:
    fig = plt.figure(figsize=(12, 4.8))
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.plot(times, states[:, 0], color='#be123c', lw=1.0)
    ax1.set_title('x(t) para r = 25')
    ax1.set_xlabel('Tiempo')
    ax1.set_ylabel('x(t)')
    ax1.grid(alpha=0.25)

    ax2 = fig.add_subplot(1, 2, 2, projection='3d')
    ax2.plot(states[:, 0], states[:, 1], states[:, 2], color='#7e22ce', lw=0.5)
    ax2.set_title('Espacio de fases para r = 25')
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    ax2.set_zlabel('z')
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def make_sensitivity_plot(times: np.ndarray, states_a: np.ndarray, states_b: np.ndarray, path: Path) -> None:
    diff = np.linalg.norm(states_a - states_b, axis=1)
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    axes[0].plot(times, states_a[:, 0], label='CI A: (1.0, 1.0, 1.0)', color='#1d4ed8', lw=1.0)
    axes[0].plot(times, states_b[:, 0], label='CI B: (1.001, 1.0, 1.0)', color='#dc2626', lw=1.0, alpha=0.8)
    axes[0].set_title('Sensibilidad a condiciones iniciales para r = 30')
    axes[0].set_ylabel('x(t)')
    axes[0].legend(frameon=False)
    axes[0].grid(alpha=0.25)

    axes[1].plot(times, diff, color='#111827', lw=1.0)
    axes[1].set_yscale('log')
    axes[1].set_xlabel('Tiempo')
    axes[1].set_ylabel('||Delta y||')
    axes[1].grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def make_r30_phase_plot(states: np.ndarray, path: Path) -> None:
    fig = plt.figure(figsize=(7.5, 5.8))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(states[:, 0], states[:, 1], states[:, 2], color='#0f172a', lw=0.5)
    ax.set_title('Espacio de fases para r = 30')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def make_r30_phase_comparison_plot(states_a: np.ndarray, states_b: np.ndarray, path: Path) -> None:
    fig = plt.figure(figsize=(12, 5.8))

    ax1 = fig.add_subplot(1, 2, 1, projection='3d')
    ax1.plot(states_a[:, 0], states_a[:, 1], states_a[:, 2], color='#1d4ed8', lw=0.5)
    ax1.set_title('CI A: (1.0, 1.0, 1.0)')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_zlabel('z')

    ax2 = fig.add_subplot(1, 2, 2, projection='3d')
    ax2.plot(states_b[:, 0], states_b[:, 1], states_b[:, 2], color='#dc2626', lw=0.5)
    ax2.set_title('CI B: (1.001, 1.0, 1.0)')
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    ax2.set_zlabel('z')

    fig.suptitle('Comparacion de trayectorias en espacio de fases para r = 30', y=0.98)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def build_summary(cases: dict[float, tuple[np.ndarray, np.ndarray]]) -> dict[str, object]:
    summary: dict[str, object] = {
        'parameters': {
            'sigma': SIGMA,
            'beta': BETA,
            'dt': DT,
            't_stop': T_STOP,
            'initial_condition': INITIAL.tolist(),
            'perturbed_initial_condition': INITIAL_PERTURBED.tolist(),
        },
        'cases': {},
    }
    for r, (times, states) in cases.items():
        x = states[:, 0]
        y = states[:, 1]
        z = states[:, 2]
        summary['cases'][str(r)] = {
            'final_state': states[-1].tolist(),
            'x_mean_last_10_time_units': float(x[times >= (T_STOP - 10.0)].mean()),
            'x_std_last_10_time_units': float(x[times >= (T_STOP - 10.0)].std()),
            'tail_std_x': settling_window_std(x),
            'tail_std_y': settling_window_std(y),
            'tail_std_z': settling_window_std(z),
            'min_distance_nonzero_fixed_point': min_distance_to_nonzero_fixed_point(states, r),
            'x_min': float(x.min()),
            'x_max': float(x.max()),
        }

    times_30, states_30 = cases[30.0]
    _, states_30_b = integrate_lorenz(30.0, INITIAL_PERTURBED)
    summary['sensitivity_r30'] = {
        'divergence_time_norm_gt_5': divergence_time(states_30, states_30_b, times_30, threshold=5.0),
        'divergence_time_norm_gt_10': divergence_time(states_30, states_30_b, times_30, threshold=10.0),
        'final_distance': float(np.linalg.norm(states_30[-1] - states_30_b[-1])),
    }
    return summary


def build_notebook() -> nbformat.NotebookNode:
    nb = nbformat.v4.new_notebook()
    nb.metadata['kernelspec'] = {
        'display_name': 'Python 3',
        'language': 'python',
        'name': 'python3',
    }
    nb.metadata['language_info'] = {'name': 'python', 'version': '3'}
    cells = []

    cells.append(
        nbformat.v4.new_markdown_cell(
            '# Guia practica: sistema de Lorenz con Runge-Kutta de cuarto orden\n\n'
            'Este notebook reproduce la integracion del sistema de Lorenz con un esquema `RK4` '
            'para estudiar equilibrio, transicion al caos y sensibilidad a condiciones iniciales.'
        )
    )

    cells.append(
        nbformat.v4.new_code_cell(
            'import json\n'
            'from pathlib import Path\n\n'
            'import matplotlib.pyplot as plt\n'
            'import numpy as np\n\n'
            "BASE_DIR = Path.cwd().resolve().parents[0]\n"
            "SUMMARY_PATH = BASE_DIR / 'content' / 'files_lorenz_rk4' / 'resumen_lorenz_rk4.json'\n"
            "with SUMMARY_PATH.open('r', encoding='utf-8') as f:\n"
            '    summary = json.load(f)\n\n'
            "summary['parameters']"
        )
    )

    cells.append(
        nbformat.v4.new_markdown_cell(
            '## Sistema y esquema numerico\n\n'
            'Se resolvio el sistema\n\n'
            '$$\\dot{x} = \\sigma (y-x), \\quad \\dot{y} = x(r-z) - y, \\quad \\dot{z} = xy - bz$$\n\n'
            'con `sigma = 10`, `b = 8/3`, `dt = 0.005`, `t in [0, 50]` y condicion inicial '
            '`(x0, y0, z0) = (1, 1, 1)`.'
        )
    )

    for title, image_name in [
        ('## Caso r = 2', 'lorenz_r2_xt.png'),
        ('## Comparacion r = 10 y r = 24: evolucion temporal', 'lorenz_r10_r24_xt.png'),
        ('## Comparacion r = 10 y r = 24: espacio de fases', 'lorenz_r10_r24_fases.png'),
        ('## Caso r = 25', 'lorenz_r25_xt_fases.png'),
        ('## Caso r = 30: espacio de fases', 'lorenz_r30_fases.png'),
        ('## Caso r = 30: comparacion de trayectorias en fases', 'lorenz_r30_fases_comparacion.png'),
        ('## Caso r = 30: sensibilidad a condiciones iniciales', 'lorenz_r30_sensibilidad.png'),
    ]:
        cells.append(nbformat.v4.new_markdown_cell(title))
        cells.append(
            nbformat.v4.new_code_cell(
                "image_path = BASE_DIR / 'content' / 'files_lorenz_rk4' / '" + image_name + "'\n"
                'img = plt.imread(image_path)\n'
                'plt.figure(figsize=(11, 5.5))\n'
                'plt.imshow(img)\n'
                'plt.axis(\'off\')\n'
                'plt.show()'
            )
        )

    cells.append(
        nbformat.v4.new_markdown_cell('## Resumen numerico')
    )
    cells.append(
        nbformat.v4.new_code_cell(
            'rows = []\n'
            "for key, values in summary['cases'].items():\n"
            '    rows.append({\n'
            "        'r': key,\n"
            "        'x_std_ultimos_10': round(values['x_std_last_10_time_units'], 4),\n"
            "        'dist_min_equilibrio_no_trivial': None if values['min_distance_nonzero_fixed_point'] is None else round(values['min_distance_nonzero_fixed_point'], 4),\n"
            "        'x_min': round(values['x_min'], 4),\n"
            "        'x_max': round(values['x_max'], 4),\n"
            '    })\n'
            'rows'
        )
    )

    cells.append(
        nbformat.v4.new_code_cell(
            "summary['sensitivity_r30']"
        )
    )

    nb.cells = cells
    return nb


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    cases = {
        2.0: integrate_lorenz(2.0, INITIAL),
        10.0: integrate_lorenz(10.0, INITIAL),
        24.0: integrate_lorenz(24.0, INITIAL),
        25.0: integrate_lorenz(25.0, INITIAL),
        30.0: integrate_lorenz(30.0, INITIAL),
    }
    times_30, states_30 = cases[30.0]
    _, states_30_b = integrate_lorenz(30.0, INITIAL_PERTURBED)

    make_time_plot(*cases[2.0], 2.0, FIG_DIR / 'lorenz_r2_xt.png')
    make_comparison_time_plot(cases, FIG_DIR / 'lorenz_r10_r24_xt.png')
    make_comparison_phase_plot(cases, FIG_DIR / 'lorenz_r10_r24_fases.png')
    make_r25_combo_plot(*cases[25.0], FIG_DIR / 'lorenz_r25_xt_fases.png')
    make_r30_phase_plot(states_30, FIG_DIR / 'lorenz_r30_fases.png')
    make_r30_phase_comparison_plot(states_30, states_30_b, FIG_DIR / 'lorenz_r30_fases_comparacion.png')
    make_sensitivity_plot(times_30, states_30, states_30_b, FIG_DIR / 'lorenz_r30_sensibilidad.png')
    make_phase_plot(cases[10.0][1], 10.0, FIG_DIR / 'lorenz_r10_fases.png')
    make_phase_plot(cases[24.0][1], 24.0, FIG_DIR / 'lorenz_r24_fases.png')
    make_phase_plot(cases[25.0][1], 25.0, FIG_DIR / 'lorenz_r25_fases.png')

    summary = build_summary(cases)
    with SUMMARY_PATH.open('w', encoding='utf-8') as f:
        json.dump(summary, f, indent=2)

    nb = build_notebook()
    with NOTEBOOK_PATH.open('w', encoding='utf-8') as f:
        nbformat.write(nb, f)

    client = NotebookClient(nb, timeout=120, kernel_name='python3')
    client.execute(cwd=str(NOTEBOOK_PATH.parent))
    with NOTEBOOK_PATH.open('w', encoding='utf-8') as f:
        nbformat.write(nb, f)

    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
