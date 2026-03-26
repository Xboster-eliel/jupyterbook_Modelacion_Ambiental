#!/usr/bin/env bash
set -euo pipefail

BOOK_DIR="/home/gentek-g3-esp/jupyterbook/Temas/modelacion ambi/book_test"
VENV_BIN="${BOOK_DIR}/.venv/bin"

cd "${BOOK_DIR}"

if [ -x "${VENV_BIN}/python" ]; then
  "${VENV_BIN}/python" -m jupyter_book start .
elif command -v jupyter-book >/dev/null 2>&1; then
  jupyter-book start .
elif command -v python3 >/dev/null 2>&1; then
  python3 -m jupyter_book start .
else
  echo "No se encontro una instalacion util de jupyter-book." >&2
  echo "Active un entorno con jupyter-book o instale dependencias con: pip install -r ../requirements.txt" >&2
  exit 1
fi
