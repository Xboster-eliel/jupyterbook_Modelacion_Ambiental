from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
from matplotlib import animation
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
import nbformat
from nbclient import NotebookClient
import numpy as np

PR = 10.0
BETA = 8.0 / 3.0
DT = 0.005
T_STOP = 50.0
INITIAL = np.array([0.0, 0.5, 0.5], dtype=float)
INITIAL_PERTURBED = np.array([0.0, 0.5, 0.50001], dtype=float)
BASE_DIR = Path(__file__).resolve().parents[1]
CONTENT_DIR = BASE_DIR / 'content'
FIG_DIR = CONTENT_DIR / 'files_lorenz_rk4'
NOTEBOOK_PATH = BASE_DIR / 'notebooks' / 'Guia_practica_Lorenz_RK4.ipynb'
SUMMARY_PATH = FIG_DIR / 'resumen_lorenz_rk4.json'


def lorenz_rhs(state: np.ndarray, r: float, pr: float = PR, beta: float = BETA) -> np.ndarray:
    w, t1, t2 = state
    return np.array(
        [
            pr * (t1 - w),
            -w * t2 + r * w - t1,
            w * t1 - beta * t2,
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


def divergence_time(states_a: np.ndarray, states_b: np.ndarray, times: np.ndarray, threshold: float = 1.0) -> float | None:
    distances = np.linalg.norm(states_a - states_b, axis=1)
    hits = np.where(distances >= threshold)[0]
    if hits.size == 0:
        return None
    return float(times[hits[0]])


def tail_std(series: np.ndarray, window: int = 2000) -> float:
    return float(np.std(series[-window:]))


def make_t1_t2_plot(times: np.ndarray, states: np.ndarray, r: float, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(9.5, 4.8))
    ax.plot(times, states[:, 1], label='T1', color='#0f766e', lw=1.1)
    ax.plot(times, states[:, 2], label='T2', color='#b45309', lw=1.1)
    ax.set_title(f'Evolucion temporal de T1 y T2 para r = {r}')
    ax.set_xlabel('Tiempo')
    ax.set_ylabel('Amplitud')
    ax.legend(frameon=False)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def make_phase_plot(states: np.ndarray, r: float, path: Path) -> None:
    fig = plt.figure(figsize=(7.2, 5.6))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(states[:, 0], states[:, 1], states[:, 2], color='#7c3aed', lw=0.55)
    ax.set_title(f'Trayectoria en espacio de fases para r = {r}')
    ax.set_xlabel('W')
    ax.set_ylabel('T1')
    ax.set_zlabel('T2')
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def make_matlab_style_collage(times: np.ndarray, states: np.ndarray, r: float, path: Path) -> None:
    fig = plt.figure(figsize=(18.0, 10.0))
    fig.patch.set_facecolor('#d9d9d9')
    gs = fig.add_gridspec(
        3,
        3,
        width_ratios=[1.08, 1.08, 1.45],
        left=0.045,
        right=0.985,
        bottom=0.06,
        top=0.90,
        wspace=0.18,
        hspace=0.42,
    )

    temporal_specs = [
        (0, states[:, 0], 'Temporal W', 'Estado W', '#1d4ed8'),
        (1, states[:, 1], 'Temporal T1', 'Estado T1', '#1d4ed8'),
        (2, states[:, 2], 'Temporal T2', 'Estado T2', '#1d4ed8'),
    ]
    for row, series, title, ylabel, color in temporal_specs:
        ax = fig.add_subplot(gs[row, 0])
        ax.plot(times, series, color=color, lw=0.45)
        ax.set_title(title, fontsize=10, pad=8)
        ax.set_xlabel('Tiempo')
        ax.set_ylabel(ylabel)
        ax.grid(True, linestyle='--', linewidth=0.6, color='k', alpha=0.55)

    phase_specs = [
        (0, states[:, 0], states[:, 1], 'Plano de fase W-T1', 'Estado W', 'Estado T1'),
        (1, states[:, 0], states[:, 2], 'Plano de fase W-T2', 'Estado W', 'Estado T2'),
        (2, states[:, 1], states[:, 2], 'Plano de fase T1-T2', 'Estado T1', 'Estado T2'),
    ]
    phase_colors = ['r', 'm', 'b']
    for (row, xdata, ydata, title, xlabel, ylabel), color in zip(phase_specs, phase_colors):
        ax = fig.add_subplot(gs[row, 1])
        ax.plot(xdata, ydata, color=color, alpha=0.7, linewidth=0.3)
        ax.set_title(title, fontsize=10, pad=8)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid(True, linestyle='--', linewidth=0.6, color='k', alpha=0.55)

    ax3d = fig.add_subplot(gs[:, 2], projection='3d')
    ax3d.plot(states[:, 0], states[:, 1], states[:, 2], color='#1d4ed8', lw=0.45)
    ax3d.set_title('Fase 3D W-T1-T2', fontsize=10, pad=16)
    ax3d.set_xlabel('Estado W', labelpad=8)
    ax3d.set_ylabel('Estado T1', labelpad=8)
    ax3d.set_zlabel('Estado T2', labelpad=8)
    ax3d.view_init(elev=26, azim=-58)
    ax3d.set_box_aspect((1.15, 1.0, 1.35))
    ax3d.grid(True)

    fig.suptitle(f'Comportamiento del sistema para r = {r}', fontsize=14, y=0.965)
    fig.savefig(path, dpi=180, facecolor=fig.get_facecolor(), bbox_inches='tight')
    plt.close(fig)


def make_t1_t2_comparison_plot(cases: dict[float, tuple[np.ndarray, np.ndarray]], r_values: list[float], path: Path) -> None:
    fig, axes = plt.subplots(len(r_values), 1, figsize=(10, 4.2 * len(r_values)), sharex=True)
    if len(r_values) == 1:
        axes = [axes]
    colors = {'T1': '#0f766e', 'T2': '#b45309'}
    for ax, r in zip(axes, r_values):
        times, states = cases[r]
        ax.plot(times, states[:, 1], color=colors['T1'], lw=1.0, label='T1')
        ax.plot(times, states[:, 2], color=colors['T2'], lw=1.0, label='T2')
        ax.set_title(f'T1(t) y T2(t) para r = {int(r)}')
        ax.set_ylabel('Amplitud')
        ax.grid(alpha=0.25)
        ax.legend(frameon=False, ncol=2, loc='upper right')
    axes[-1].set_xlabel('Tiempo')
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def make_phase_comparison_plot(cases: dict[float, tuple[np.ndarray, np.ndarray]], r_values: list[float], path: Path) -> None:
    fig = plt.figure(figsize=(6.2 * len(r_values), 5.5))
    colors = ['#0f766e', '#b45309', '#7e22ce']
    for i, r in enumerate(r_values, start=1):
        ax = fig.add_subplot(1, len(r_values), i, projection='3d')
        states = cases[r][1]
        ax.plot(states[:, 0], states[:, 1], states[:, 2], color=colors[i - 1], lw=0.5)
        ax.set_title(f'Fases para r = {int(r)}')
        ax.set_xlabel('W')
        ax.set_ylabel('T1')
        ax.set_zlabel('T2')
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def make_t1_sensitivity_plot(times: np.ndarray, states_a: np.ndarray, states_b: np.ndarray, path: Path) -> None:
    diff = np.linalg.norm(states_a - states_b, axis=1)
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    axes[0].plot(times, states_a[:, 1], label='CI A: (0, 0.5, 0.5)', color='#1d4ed8', lw=1.0)
    axes[0].plot(times, states_b[:, 1], label='CI B: (0, 0.5, 0.50001)', color='#dc2626', lw=1.0, alpha=0.85)
    axes[0].set_title('Evolucion temporal de T1 para r = 30')
    axes[0].set_ylabel('T1')
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


def make_r30_phase_comparison_plot(states_a: np.ndarray, states_b: np.ndarray, path: Path) -> None:
    fig = plt.figure(figsize=(12, 5.8))

    ax1 = fig.add_subplot(1, 2, 1, projection='3d')
    ax1.plot(states_a[:, 0], states_a[:, 1], states_a[:, 2], color='#1d4ed8', lw=0.5)
    ax1.set_title('CI A: (0, 0.5, 0.5)')
    ax1.set_xlabel('W')
    ax1.set_ylabel('T1')
    ax1.set_zlabel('T2')

    ax2 = fig.add_subplot(1, 2, 2, projection='3d')
    ax2.plot(states_b[:, 0], states_b[:, 1], states_b[:, 2], color='#dc2626', lw=0.5)
    ax2.set_title('CI B: (0, 0.5, 0.50001)')
    ax2.set_xlabel('W')
    ax2.set_ylabel('T1')
    ax2.set_zlabel('T2')

    fig.suptitle('Comparacion de trayectorias en espacio de fases para r = 30', y=0.98)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def make_r30_phase_animation(
    states_a: np.ndarray,
    states_b: np.ndarray,
    path: Path,
    step: int = 16,
    tail: int = 180,
) -> None:
    fig = plt.figure(figsize=(8.6, 6.4), facecolor='white')
    ax = fig.add_subplot(111, projection='3d')

    all_states = np.vstack([states_a, states_b])
    w_pad = 2.0
    t1_pad = 2.0
    t2_pad = 2.0
    ax.set_xlim(all_states[:, 0].min() - w_pad, all_states[:, 0].max() + w_pad)
    ax.set_ylim(all_states[:, 1].min() - t1_pad, all_states[:, 1].max() + t1_pad)
    ax.set_zlim(max(0.0, all_states[:, 2].min() - t2_pad), all_states[:, 2].max() + t2_pad)
    ax.set_xlabel('Estado W', labelpad=10)
    ax.set_ylabel('Estado T1', labelpad=10)
    ax.set_zlabel('Estado T2', labelpad=10)
    ax.set_title('Animacion del atractor para r = 30', pad=14)
    ax.view_init(elev=28, azim=-52)
    ax.set_box_aspect((1.15, 1.0, 1.25))
    ax.grid(True)

    color_a = '#1d4ed8'
    color_b = '#dc2626'
    line_a, = ax.plot([], [], [], '-', color=color_a, lw=1.0, alpha=0.85, label='CI A: (0, 0.5, 0.5)')
    line_b, = ax.plot([], [], [], '-', color=color_b, lw=1.0, alpha=0.85, label='CI B: (0, 0.5, 0.50001)')
    pt_a, = ax.plot([], [], [], 'o', color=color_a, ms=5)
    pt_b, = ax.plot([], [], [], 'o', color=color_b, ms=5)
    ax.legend(loc='upper left', frameon=True, fontsize=8)

    frame_ids = list(range(2, states_a.shape[0], step))

    def init() -> list:
        for artist in (line_a, line_b, pt_a, pt_b):
            artist.set_data([], [])
            artist.set_3d_properties([])
        return [line_a, line_b, pt_a, pt_b]

    def animate(frame_index: int) -> list:
        end = frame_ids[frame_index]
        start = max(0, end - tail)

        wa, t1a, t2a = states_a[start:end].T
        wb, t1b, t2b = states_b[start:end].T

        line_a.set_data(wa, t1a)
        line_a.set_3d_properties(t2a)
        line_b.set_data(wb, t1b)
        line_b.set_3d_properties(t2b)

        pt_a.set_data([wa[-1]], [t1a[-1]])
        pt_a.set_3d_properties([t2a[-1]])
        pt_b.set_data([wb[-1]], [t1b[-1]])
        pt_b.set_3d_properties([t2b[-1]])

        ax.view_init(elev=28, azim=-52 + 0.12 * frame_index)
        return [line_a, line_b, pt_a, pt_b]

    anim = animation.FuncAnimation(
        fig,
        animate,
        init_func=init,
        frames=len(frame_ids),
        interval=120,
        blit=False,
    )
    writer = animation.PillowWriter(fps=8)
    anim.save(path, writer=writer, dpi=75)
    plt.close(fig)


def make_r30_multi_trajectory_animation(path: Path) -> None:
    n_trajectories = 5
    dt_anim = 0.01
    num_steps = 3000
    times = np.linspace(0.0, dt_anim * num_steps, num_steps + 1)

    initials = np.array(
        [[0.0, 0.5 + 0.06 * i, 0.5 + 0.012 * i] for i in range(n_trajectories)],
        dtype=float,
    )
    trajectories = np.asarray([integrate_lorenz(30.0, y0, dt=dt_anim, t_stop=times[-1])[1] for y0 in initials])

    fig = plt.figure(figsize=(8.8, 6.5), facecolor='white')
    ax = fig.add_subplot(111, projection='3d')

    all_states = trajectories.reshape(-1, 3)
    ax.set_xlim(all_states[:, 0].min() - 2.0, all_states[:, 0].max() + 2.0)
    ax.set_ylim(all_states[:, 1].min() - 2.0, all_states[:, 1].max() + 2.0)
    ax.set_zlim(max(0.0, all_states[:, 2].min() - 2.0), all_states[:, 2].max() + 2.0)
    ax.set_xlabel('Estado W', labelpad=10)
    ax.set_ylabel('Estado T1', labelpad=10)
    ax.set_zlabel('Estado T2', labelpad=10)
    ax.set_title('Animacion de cinco trayectorias cercanas para r = 30', pad=14)
    ax.view_init(elev=28, azim=-55)
    ax.set_box_aspect((1.15, 1.0, 1.25))
    ax.grid(True)

    colors = plt.cm.tab10(np.linspace(0, 1, n_trajectories))
    lines = []
    points = []
    for color, y0 in zip(colors, initials):
        line, = ax.plot([], [], [], '-', color=color, lw=0.8, alpha=0.8, label=f'y_0=({y0[0]:.2f}, {y0[1]:.2f}, {y0[2]:.3f})')
        point, = ax.plot([], [], [], 'o', color=color, ms=4.5)
        lines.append(line)
        points.append(point)
    ax.legend(loc='upper left', frameon=True, fontsize=7)

    frame_ids = list(range(2, trajectories.shape[1], 14))
    tail = 180

    def init() -> list:
        for artist in [*lines, *points]:
            artist.set_data([], [])
            artist.set_3d_properties([])
        return [*lines, *points]

    def animate(frame_index: int) -> list:
        end = frame_ids[frame_index]
        start = max(0, end - tail)
        for line, point, traj in zip(lines, points, trajectories):
            w, t1, t2 = traj[start:end].T
            line.set_data(w, t1)
            line.set_3d_properties(t2)
            point.set_data([w[-1]], [t1[-1]])
            point.set_3d_properties([t2[-1]])
        ax.view_init(elev=28, azim=-55 + 0.14 * frame_index)
        return [*lines, *points]

    anim = animation.FuncAnimation(
        fig,
        animate,
        init_func=init,
        frames=len(frame_ids),
        interval=120,
        blit=False,
    )
    writer = animation.PillowWriter(fps=8)
    anim.save(path, writer=writer, dpi=75)
    plt.close(fig)


def build_summary(cases: dict[float, tuple[np.ndarray, np.ndarray]]) -> dict[str, object]:
    summary: dict[str, object] = {
        'parameters': {
            'Pr': PR,
            'b': BETA,
            'dt': DT,
            't_stop': T_STOP,
            'initial_condition': INITIAL.tolist(),
            'perturbed_initial_condition': INITIAL_PERTURBED.tolist(),
        },
        'cases': {},
    }
    for r, (times, states) in cases.items():
        w = states[:, 0]
        t1 = states[:, 1]
        t2 = states[:, 2]
        summary['cases'][str(r)] = {
            'final_state': states[-1].tolist(),
            't1_mean_last_10_time_units': float(t1[times >= (T_STOP - 10.0)].mean()),
            't1_std_last_10_time_units': float(t1[times >= (T_STOP - 10.0)].std()),
            't2_std_last_10_time_units': float(t2[times >= (T_STOP - 10.0)].std()),
            'tail_std_w': tail_std(w),
            'tail_std_t1': tail_std(t1),
            'tail_std_t2': tail_std(t2),
            'min_distance_nonzero_fixed_point': min_distance_to_nonzero_fixed_point(states, r),
            'w_min': float(w.min()),
            'w_max': float(w.max()),
            't1_min': float(t1.min()),
            't1_max': float(t1.max()),
            't2_min': float(t2.min()),
            't2_max': float(t2.max()),
        }

    times_30, states_30 = cases[30.0]
    _, states_30_b = integrate_lorenz(30.0, INITIAL_PERTURBED)
    summary['sensitivity_r30'] = {
        'divergence_time_norm_gt_1': divergence_time(states_30, states_30_b, times_30, threshold=1.0),
        'divergence_time_norm_gt_5': divergence_time(states_30, states_30_b, times_30, threshold=5.0),
        'final_distance': float(np.linalg.norm(states_30[-1] - states_30_b[-1])),
    }
    return summary


def image_cell(image_name: str) -> nbformat.NotebookNode:
    return nbformat.v4.new_code_cell(
        "image_path = BASE_DIR / 'content' / 'files_lorenz_rk4' / '" + image_name + "'\n"
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
            'Este notebook resuelve el sistema del inciso con la notacion `W`, `T1` y `T2`, '
            'usando `y_0 = (0, 0.5, 0.5)`, `Pr = 10`, `b = 8/3`, `dt = 0.005` y `0 <= t <= 50`. '
            'El objetivo no es solo reproducir figuras, sino interpretar criticamente la transicion '
            'desde un equilibrio estable hacia un comportamiento caotico.'
        )
    )

    cells.append(
        nbformat.v4.new_code_cell(
            'import json\n'
            'from pathlib import Path\n\n'
            'import matplotlib.pyplot as plt\n\n'
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
            '$$\\frac{dW}{dt} = Pr (T_1 - W), \\quad '
            '\\frac{dT_1}{dt} = -W T_2 + rW - T_1, \\quad '
            '\\frac{dT_2}{dt} = W T_1 - bT_2$$\n\n'
            'con `Pr = 10`, `b = 8/3`, `dt = 0.005` y condicion inicial base `y_0 = (0, 0.5, 0.5)`. '
            'La implementacion numerica sigue la estructura clasica de `RK4`; la lectura del problema se centra '
            'en comparar la evolucion temporal de `T1(t)` y `T2(t)` con la trayectoria en el espacio de fases `W-T1-T2`.'
        )
    )

    cells.append(nbformat.v4.new_markdown_cell('## Caso `r = 2`: evolucion temporal de `T1(t)` y `T2(t)`'))
    cells.append(image_cell('lorenz_r2_t1_t2.png'))
    cells.append(
        nbformat.v4.new_markdown_cell(
            'La lectura clave de este caso es la amortiguacion. `T1` y `T2` pierden variabilidad y convergen '
            'hacia valores casi constantes, lo que indica un equilibrio estable. Este caso sirve como referencia '
            'para distinguir despues los regimens transicionales y caoticos.'
        )
    )

    cells.append(nbformat.v4.new_markdown_cell('## Casos `r = 10` y `r = 24`: evolucion temporal de `T1(t)` y `T2(t)`'))
    cells.append(image_cell('lorenz_r10_r24_t1_t2.png'))
    cells.append(nbformat.v4.new_markdown_cell('## Caso `r = 10`: collage temporal-fases estilo MATLAB'))
    cells.append(image_cell('lorenz_r10_collage_matlab.png'))
    cells.append(nbformat.v4.new_markdown_cell('## Caso `r = 24`: collage temporal-fases estilo MATLAB'))
    cells.append(image_cell('lorenz_r24_collage_matlab.png'))
    cells.append(
        nbformat.v4.new_markdown_cell(
            'La comparacion entre estos dos casos es crucial. El collage, ahora fiel a la organizacion usada en MATLAB, '
            'permite leer a la vez las tres series temporales, las tres proyecciones de fase 2D y la trayectoria 3D. '
            'En `r = 10`, todas esas vistas siguen concentrandose alrededor de un equilibrio. En `r = 24`, en cambio, '
            'aparecen recorridos mas amplios y cambios entre regiones del espacio de fases, lo que marca un regimen transicional.'
        )
    )

    cells.append(nbformat.v4.new_markdown_cell('## Caso `r = 25`: evolucion temporal de `T1(t)` y `T2(t)`'))
    cells.append(image_cell('lorenz_r25_t1_t2.png'))
    cells.append(nbformat.v4.new_markdown_cell('## Caso `r = 25`: collage temporal-fases estilo MATLAB'))
    cells.append(image_cell('lorenz_r25_collage_matlab.png'))
    cells.append(
        nbformat.v4.new_markdown_cell(
            'En `r = 25` la interpretacion es mas nítida que en `r = 24`: el sistema ya no solo parece irregular, '
            'sino que muestra una dinamica persistentemente no periodica dentro de una estructura atractora acotada. '
            'La comparacion entre ambos casos ayuda a evitar una lectura apresurada de `r = 24` como si fuera un equilibrio muy lento.'
        )
    )

    cells.append(nbformat.v4.new_markdown_cell('## Caso `r = 30`: collage temporal-fases estilo MATLAB'))
    cells.append(image_cell('lorenz_r30_collage_matlab.png'))
    cells.append(nbformat.v4.new_markdown_cell('## Caso `r = 30`: animacion de trayectorias en el espacio de fases'))
    cells.append(markdown_image_cell('../content/files_lorenz_rk4/lorenz_r30_fases_animacion.gif', 'Animacion del atractor para r = 30'))
    cells.append(nbformat.v4.new_markdown_cell('## Caso `r = 30`: evolucion temporal de `T1(t)` y sensibilidad a condiciones iniciales'))
    cells.append(image_cell('lorenz_r30_t1_sensibilidad.png'))
    cells.append(
        nbformat.v4.new_markdown_cell(
            'La lectura conjunta del collage, de la primera animacion y de la serie temporal de `T1` permite fijar la '
            'conclusion principal del caso `r = 30`. Las dos trayectorias parten de estados casi indistinguibles, '
            'permanecen dentro del mismo atractor y, sin embargo, pierden rapidamente la coincidencia punto a punto. '
            'Esto significa que el sistema conserva una estructura geometrica global robusta, pero pierde predictibilidad puntual '
            'en horizontes de tiempo finitos, que es precisamente la firma del caos determinista.'
        )
    )
    cells.append(nbformat.v4.new_markdown_cell('## Caso `r = 30`: animacion de cinco trayectorias cercanas'))
    cells.append(markdown_image_cell('../content/files_lorenz_rk4/lorenz_r30_multitrayectorias_animacion.gif', 'Animacion de cinco trayectorias cercanas para r = 30'))
    cells.append(
        nbformat.v4.new_markdown_cell(
            'La segunda animacion amplia la interpretacion anterior desde una sola pareja de trayectorias hacia un pequeño conjunto '
            'de orbitas vecinas. Su interes analitico es mostrar que la sensibilidad a condiciones iniciales no depende de una eleccion '
            'especialmente desafortunada de dos estados iniciales, sino que es una propiedad estructural del regimen caotico para `r = 30`. '
            'Las cinco trayectorias quedan confinadas dentro de la misma geometria global del atractor, pero recorren con rapidez regiones '
            'distintas del espacio de fases y distribuyen de forma desigual sus tiempos de permanencia en cada lobo. Esta observacion es clave '
            'porque anticipa una idea mas profunda: en sistemas caoticos, la trayectoria puntual deja de ser reproducible a largo plazo, '
            'pero la estructura estadistica del atractor sigue siendo interpretable y comparable mediante observables resumidos.'
        )
    )

    cells.append(nbformat.v4.new_markdown_cell('## Resumen numerico'))
    cells.append(
        nbformat.v4.new_code_cell(
            'rows = []\n'
            "for key, values in summary['cases'].items():\n"
            '    rows.append({\n'
            "        'r': key,\n"
            "        'std_T1_ultimos_10': round(values['t1_std_last_10_time_units'], 4),\n"
            "        'std_T2_ultimos_10': round(values['t2_std_last_10_time_units'], 4),\n"
            "        'dist_min_equilibrio_no_trivial': None if values['min_distance_nonzero_fixed_point'] is None else round(values['min_distance_nonzero_fixed_point'], 4),\n"
            "        'T1_min': round(values['t1_min'], 4),\n"
            "        'T1_max': round(values['t1_max'], 4),\n"
            '    })\n'
            'rows'
        )
    )
    cells.append(nbformat.v4.new_code_cell("summary['sensitivity_r30']"))
    cells.append(
        nbformat.v4.new_markdown_cell(
            '## Cierre interpretativo\n\n'
            'El patron general del capitulo es claro: al aumentar `r`, el sistema pasa de amortiguar perturbaciones '\
            'a sostener oscilaciones persistentes y, finalmente, a mostrar sensibilidad extrema a condiciones iniciales. '\
            'La implementacion de `RK4` permite capturar bien esa transicion, pero el aprendizaje principal proviene de '\
            'comparar criticamente las salidas temporales y el espacio de fases en conjunto.'
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

    make_t1_t2_plot(*cases[2.0], 2.0, FIG_DIR / 'lorenz_r2_t1_t2.png')
    make_t1_t2_comparison_plot(cases, [10.0, 24.0], FIG_DIR / 'lorenz_r10_r24_t1_t2.png')
    make_matlab_style_collage(*cases[10.0], 10.0, FIG_DIR / 'lorenz_r10_collage_matlab.png')
    make_matlab_style_collage(*cases[24.0], 24.0, FIG_DIR / 'lorenz_r24_collage_matlab.png')
    make_t1_t2_plot(*cases[25.0], 25.0, FIG_DIR / 'lorenz_r25_t1_t2.png')
    make_matlab_style_collage(*cases[25.0], 25.0, FIG_DIR / 'lorenz_r25_collage_matlab.png')
    make_matlab_style_collage(*cases[30.0], 30.0, FIG_DIR / 'lorenz_r30_collage_matlab.png')
    make_r30_phase_animation(states_30, states_30_b, FIG_DIR / 'lorenz_r30_fases_animacion.gif')
    make_r30_multi_trajectory_animation(FIG_DIR / 'lorenz_r30_multitrayectorias_animacion.gif')
    make_t1_sensitivity_plot(times_30, states_30, states_30_b, FIG_DIR / 'lorenz_r30_t1_sensibilidad.png')

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
