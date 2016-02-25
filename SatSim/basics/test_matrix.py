import unittest
from vector_types import Matrix


class TestMatrix(unittest.TestCase):
    def test_transpose(self):
        m = Matrix([[1, 2, 3],
                    [4, 5, 6],
                    [7, 8, 9]])
        m_transpose = m.T()
        self.assertTrue(m_transpose == Matrix([[1, 4, 7], [2, 5, 8], [3, 6, 9]]))

    def test_addition(self):
        m1 = Matrix([1.5, 2.0, 2, 12.0, 11, 13.0, 23, -12, -123.0])
        m2 = Matrix([2.0, -30, 13, 2.4, -90.5, -10.4, 12.3, 12, 1])
        self.assertAlmostEqual(m1+m2, Matrix([3.5, -28, 15, 14.4, -79.5, 2.6, 35.3, 0, -122]))

    def test_subtraction(self):
        m1 = Matrix([1, 2, 3, 4, 5, 6, 7, 8, 9])
        m2 = Matrix([9, 8, 7, 6, 5, 4, 3, 2, 1])
        self.assertAlmostEqual(m1-m2, Matrix([-8, -6, -4, -2, 0, 2, 4, 6, 8]))

    def test_matrix_multiplication(self):
        m1 = Matrix([1, 0, 0, 0, 1, 0, 0, 0, 1])
        m2 = Matrix([2, 5, -4, 2, -9, 3, -4.5, -11, 3])
        self.assertAlmostEqual(m1*m2, m2)

    def test_scalar_multiplication(self):
        m = Matrix([-1.0, -2.0, -3.0, -4, -5, -6, -7, -8, -9])
        self.assertEqual(m*(-1), Matrix([1, 2, 3, 4, 5, 6, 7, 8, 9]))

    def test_scalar_division(self):
        m = Matrix([2, 4, 6, 8, 10, 12, 14, 16, 18])
        self.assertAlmostEqual(m/2, Matrix([1, 2, 3, 4, 5, 6, 7, 8, 9]))


if __name__ == "__main__":
    unittest.main()
