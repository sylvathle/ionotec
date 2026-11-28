from pathlib import Path
from typing import Iterable, Union
import shutil
import gzip
import subprocess

def extract_and_remove_gz(gz_path: Path) -> Path:
    extracted_path = gz_path.with_suffix("")

    with gzip.open(gz_path, "rb") as f_in, open(extracted_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)

    gz_path.unlink()

    return extracted_path

def remote_exists(host, path):
    r = subprocess.run(["ssh", host, f"test -f {path}"], capture_output=True)
    if r.returncode == 255:
        raise ConnectionError(r.stderr.decode().strip())
    return r.returncode == 0

def delete_files(paths: Iterable[Union[str, "Path"]], missing_ok: bool = True) -> None:
    """Delete each file in `paths`. Ignores missing files if missing_ok=True."""
    for p in paths:
        p = Path(p)
        try:
            p.unlink(missing_ok=missing_ok)
        except IsADirectoryError:
            print(f"Skipping directory: {p}")
        except PermissionError as e:
            print(f"Could not delete {p}: {e}")
