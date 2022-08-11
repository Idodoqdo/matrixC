#include "s21_matrix.h"

matrix_t s21_create_matrix(int rows, int columns) {
    matrix_t A;
    if (rows < 0 || columns < 0) {
        A.matrix_type = INCORRECT_MATRIX;
        A.columns = 0;
        A.rows = 0;
    } else {
    A.matrix = (double**)calloc(rows, sizeof(double*));
    A.rows = rows;
    A.columns = columns;
    for (int i = 0; i < rows; i++) {
        A.matrix[i] = (double*)calloc(columns, sizeof(double));
    }
    s21_tip(&A);
}
return A;
}

void s21_tip(matrix_t *A) {
    int edmatr = 0;
    int zeromatr = 0;
    for (int z = 0; z < A->rows; z++) {
        for (int x = 0; x < A->columns; x++) {
            zeromatr += fabs(A->matrix[z][x]);
            if (z == x && A->matrix[z][x] == 1) {
                edmatr += A->matrix[z][x];
            }
        }
    }
    if (zeromatr == 0) {
        A->matrix_type = ZERO_MATRIX;
    } else if (edmatr == A->rows && A->columns == A->rows) {
        A->matrix_type = IDENTITY_MATRIX;
    } else {
        A->matrix_type = CORRECT_MATRIX;
    }
}

void s21_remove_matrix(matrix_t *A) {
    for (int i = 0; i < A->rows; i++) {
        free(A->matrix[i]);
    }
    free(A->matrix);
    A->columns = 0;
    A->rows = 0;
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
    int rez = SUCCESS;
    if (A->rows != B->rows || A->columns != B->columns || A->matrix_type == INCORRECT_MATRIX
    || B->matrix_type == INCORRECT_MATRIX) {
        rez = FAILURE;
    } else {
        for (int i = 0; i < A->rows; i++) {
            for (int j = 0; j < A->columns; j++) {
                if (fabs(A->matrix[i][j] - B->matrix[i][j]) >= 0.000001) {
                    rez = FAILURE;
                }
            }
        }
    }
    return rez;
}

matrix_t s21_sum_matrix(matrix_t *A, matrix_t *B) {
    matrix_t C = s21_create_matrix(A->rows, A->columns);
    if (A->rows != B->rows || A->columns != B->columns || A->matrix_type == INCORRECT_MATRIX
    || B->matrix_type == INCORRECT_MATRIX) {
        C.matrix_type = INCORRECT_MATRIX;
    } else {
        for (int i = 0; i < A->rows; i++) {
            for (int j = 0; j < A->columns; j++) {
               C.matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
            }
        }
        s21_tip(&C);
    }
    return C;
}

matrix_t s21_sub_matrix(matrix_t *A, matrix_t *B) {
    matrix_t C = s21_create_matrix(A->rows, A->columns);
    if (A->rows != B->rows || A->columns != B->columns || A->matrix_type == INCORRECT_MATRIX
    || B->matrix_type == INCORRECT_MATRIX) {
        C.matrix_type = INCORRECT_MATRIX;
    } else {
        for (int i = 0; i < A->rows; i++) {
            for (int j = 0; j < A->columns; j++) {
               C.matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
            }
        }
        s21_tip(&C);
    }
    return C;
}

matrix_t s21_mult_number(matrix_t *A, double number) {
    matrix_t B = s21_create_matrix(A->rows, A->columns);
    if (A->matrix_type == INCORRECT_MATRIX) {
         B.matrix_type = INCORRECT_MATRIX;
    } else {
        for (int i = 0; i < A->rows; i++) {
            for (int j = 0; j < A->columns; j++) {
               B.matrix[i][j] = A->matrix[i][j] * number;
            }
        }
        s21_tip(&B);
    }
    return B;
}

matrix_t s21_mult_matrix(matrix_t *A, matrix_t *B) {
    matrix_t C = s21_create_matrix(A->rows, B->columns);
    if (A->matrix_type == INCORRECT_MATRIX || B->matrix_type == INCORRECT_MATRIX || A->columns != B->rows) {
        C.matrix_type = INCORRECT_MATRIX;
    } else {
       C.rows = A->rows;
       C.columns = B->columns;
       for (int i = 0; i < C.rows; i++) {
            for (int j = 0; j < C.columns; j++) {
                C.matrix[i][j] = 0;
                for (int k = 0; k < B->rows; k++) {
                    C.matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
                }
            }
        }
        s21_tip(&C);
    }
    return C;
}

matrix_t s21_transpose(matrix_t *A) {
    matrix_t B = s21_create_matrix(A->columns, A->rows);
    if (A->matrix_type == INCORRECT_MATRIX) {
        B.matrix_type = INCORRECT_MATRIX;
    } else {
        for (int i = 0; i < A->rows; i++) {
            for (int j = 0; j < A->columns; j++) {
                B.matrix[j][i] = A->matrix[i][j];
            }
        }
        s21_tip(&B);
    }
    return B;
}

matrix_t s21_minor(int i, int j, matrix_t *A) {
    matrix_t B = s21_create_matrix(A->rows - 1, A->columns - 1);
    int str = 0;
    int stb = 0;
        for (int k = 0; k < A->rows; k++) {
            for (int m = 0; m < A->columns; m++) {
                if (k != i && m != j) {
                    B.matrix[str][stb] = A->matrix[k][m];
                    stb++;
                    if (stb == A->columns - 1) str++;
                }
            }
            stb = 0;
        }
    return B;
}

double s21_determinant(matrix_t *A) {
    double rez = 0;
    if (A->matrix_type == INCORRECT_MATRIX || A->rows != A->columns) {
        rez = NAN;
    } else {
        matrix_t B;
        if (A->columns == 1) {
            rez = A->matrix[0][0];
        } else if (A->columns == 2) {
            rez = A->matrix[0][0] * A->matrix[1][1] - A->matrix[1][0] * A->matrix[0][1];
        } else if (A->columns > 2) {
        int i = 0;
        for (int j = 0; j < A->columns; j++) {
            B = s21_minor(i, j, A);
            rez += pow(-1, i+j) * A->matrix[i][j] * s21_determinant(&B);
            s21_remove_matrix(&B);
        }
        }
    }
    return rez;
}

matrix_t s21_calc_complements(matrix_t *A) {
    double rez = 0;
    matrix_t C = s21_create_matrix(A->rows, A->columns);
    if (A->matrix_type == INCORRECT_MATRIX || A->rows != A->columns) {
        C.matrix_type = INCORRECT_MATRIX;
    } else {
        matrix_t B;
        for (int i = 0; i < A->rows; i++) {
            for (int j = 0; j < A->columns; j++) {
            rez = 0;
            B = s21_minor(i, j, A);
            rez = pow(-1, i+j) * s21_determinant(&B);
            C.matrix[i][j] = rez;
            s21_remove_matrix(&B);
            }
        }
        s21_tip(&C);
    }
    return C;
}

matrix_t s21_inverse_matrix(matrix_t *A) {
    matrix_t B;
    if (A->matrix_type == INCORRECT_MATRIX || (A->rows != A->columns)) {
        B.matrix_type = INCORRECT_MATRIX;
    } else {
        double opr = s21_determinant(A);
        if (opr == 0) {
            B.matrix_type = INCORRECT_MATRIX;
        } else {
            matrix_t C = s21_calc_complements(A);
            matrix_t D = s21_transpose(&C);
            s21_remove_matrix(&C);
            B = s21_mult_number(&D, 1/opr);
            s21_remove_matrix(&D);
            s21_tip(&B);
        }
        }
    return B;
}

