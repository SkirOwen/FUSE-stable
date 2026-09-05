from __future__ import annotations

import os
import glob


def clean_files(dir_path: str, file_ext: str) -> None:
	for file in glob.glob(os.path.join(dir_path, f"*.{file_ext}")):
		os.remove(file)
