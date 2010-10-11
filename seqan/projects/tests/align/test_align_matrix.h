#ifndef TESTS_ALIGN_TEST_ALIGN_MATRIX_H_
#define TESTS_ALIGN_TEST_ALIGN_MATRIX_H_

SEQAN_DEFINE_TEST(test_align_matrix)
{
    using namespace std;
    using namespace seqan;

	// Resize the matrix.
	Matrix<double,2> matrix1;

	setLength(matrix1, 0, 2);
	setLength(matrix1, 1, 2);
	resize(matrix1);
	value(matrix1,0,0) = 3.0;
	value(matrix1,0,1) = 3.5;
	value(matrix1,1,0) = 14.5;
	value(matrix1,1,1) = -6.0;

	Matrix<double,2> matrix2;

	setLength(matrix2, 0, 2);
	setLength(matrix2, 1, 2);
	fill(matrix2,1.0);

	// Basic 2D operations, A+B,A-B,A*a,A*B, A==B

	Matrix<double,2> matrix3;

	matrix3 = matrix1 + matrix2;
	SEQAN_ASSERT_EQ(value(matrix3,0,0), value(matrix1,0,0) + value(matrix2,0,0));
	SEQAN_ASSERT_EQ(value(matrix3,0,1), value(matrix1,0,1) + value(matrix2,0,1));
	SEQAN_ASSERT_EQ(value(matrix3,1,0), value(matrix1,1,0) + value(matrix2,1,0));
	SEQAN_ASSERT_EQ(value(matrix3,1,1), value(matrix1,1,1) + value(matrix2,1,1));
	SEQAN_ASSERT_EQ(matrix1 + matrix2, matrix2 + matrix1);
	matrix3 = matrix1 - matrix2;
	SEQAN_ASSERT_EQ(value(matrix3,0,0), value(matrix1,0,0) - value(matrix2,0,0));
	SEQAN_ASSERT_EQ(value(matrix3,0,1), value(matrix1,0,1) - value(matrix2,0,1));
	SEQAN_ASSERT_EQ(value(matrix3,1,0), value(matrix1,1,0) - value(matrix2,1,0));
	SEQAN_ASSERT_EQ(value(matrix3,1,1), value(matrix1,1,1) - value(matrix2,1,1));

    SEQAN_ASSERT_EQ(matrix1-matrix2, matrix1+(matrix2*(-1.0)));
	matrix3=matrix1*matrix2;
	matrix3=matrix1*5.0;
	SEQAN_ASSERT_EQ(matrix1 * 5.0, 5.0 * matrix1);
    
    // n-dimensional matrix
	Matrix<double> matrixN;
	setDimension(matrixN,2);
	setLength(matrix2, 0, 3);
	setLength(matrix2, 1, 2);
	fill(matrixN,1.0);
}

#endif  // TESTS_ALIGN_TEST_ALIGN_MATRIX_H_
