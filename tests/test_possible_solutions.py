from fuse202.structure.possible_solutions import (
	CrystalSolutionsCalculator,
	cube_function,
	orthorhombic_function,
	return_factors,
	tetragonal_function,
)


def test_return_factors_of_twelve():
	assert return_factors(12) == {1, 2, 3, 4, 6, 12}


def test_return_factors_of_prime():
	assert return_factors(13) == {1, 13}


def test_cube_function():
	assert cube_function(2) == 16  # 2 * (2**3)


def test_tetragonal_function():
	assert tetragonal_function(3, 4) == 36  # 3**2 * 4


def test_orthorhombic_function():
	assert orthorhombic_function(2, 3, 4) == 24


def test_crystal_solutions_calculator_writes_cache_files_to_given_dir(tmp_path):
	calc = CrystalSolutionsCalculator(max_ax=4, solutions_dir=str(tmp_path))
	cubic, tetragonal, hexagonal, orthorhombic, monoclinic = calc.calculate_all_solutions()

	assert (tmp_path / "cubes.p").exists()
	assert (tmp_path / "tetragonal.p").exists()
	assert (tmp_path / "hexagonal.p").exists()
	assert (tmp_path / "orthorhombic.p").exists()
	assert (tmp_path / "monoclinic.p").exists()

	assert cubic[cube_function(2)] == [2, 2, 4]


def test_crystal_solutions_calculator_reuses_cache_on_second_call(tmp_path):
	calc = CrystalSolutionsCalculator(max_ax=4, solutions_dir=str(tmp_path))
	first = calc.calculate_all_solutions()

	cache_mtime = (tmp_path / "cubes.p").stat().st_mtime
	second = calc.calculate_all_solutions()

	assert (tmp_path / "cubes.p").stat().st_mtime == cache_mtime
	assert first == second


def test_crystal_solutions_calculator_restart_recomputes(tmp_path):
	calc = CrystalSolutionsCalculator(max_ax=3, solutions_dir=str(tmp_path))
	calc.calculate_all_solutions()
	assert list(tmp_path.glob("*.p"))

	# restart=True clears the cache directory first, then recomputes - should
	# not raise, and files should exist again afterward.
	result = calc.calculate_all_solutions(restart=True)
	assert list(tmp_path.glob("*.p"))
	assert result[0]  # cubic solutions dict is non-empty
