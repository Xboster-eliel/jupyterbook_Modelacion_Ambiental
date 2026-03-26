# jupyterbook_Modelacion_Ambiental

Book del curso de Modelacion Ambiental construido con MyST y Jupyter Book.

## Requisitos

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

El workflow de GitHub Actions construye y publica el sitio en GitHub Pages en cada push a `main`.

## Desarrollo local

```bash
./start-book.sh
```

## Build

```bash
python -m jupyter_book build --all --html --site --force
```

## Sitio publicado

`https://xboster-eliel.github.io/jupyterbook_Modelacion_Ambiental/`
