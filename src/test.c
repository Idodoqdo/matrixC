#include <check.h>
#include "s21_matrix.h"


START_TEST(s21_create_1) {
    matrix_t A;
    A = s21_create_matrix(3, 4);
    int res = 1, res1 = 1;
    for (int i = 0; i < A.rows; i++) {
        for (int j = 0; j < A.columns; j++) {
            if (A.matrix[i][j] != 0.0) {
                res = 0;
            }
        }
    }
    if (A.matrix_type != ZERO_MATRIX) {
        res1 = 0;
    }
    ck_assert_int_eq(res, 1);
    ck_assert_int_eq(res1, 1);
    ck_assert_int_eq(A.rows, 3);
    ck_assert_int_eq(A.columns, 4);
    s21_remove_matrix(&A);
} END_TEST

START_TEST(s21_create_2) {
    matrix_t A;
    A = s21_create_matrix(0, -1);
    int res = 1;
    if (A.matrix_type != INCORRECT_MATRIX) {
        res = 0;
    }
    ck_assert_int_eq(res, 1);
    ck_assert_int_eq(A.rows, 0);
    ck_assert_int_eq(A.columns, 0);
} END_TEST



START_TEST(s21_create_3) {
    matrix_t a = {0};
    a = s21_create_matrix(1, 1);
    ck_assert_int_eq(a.matrix_type, ZERO_MATRIX);
    matrix_t a1 = {0};
    a1 = s21_create_matrix(-1, -2);
    ck_assert_int_eq(a1.matrix_type, INCORRECT_MATRIX);
    s21_remove_matrix(&a);
} END_TEST

START_TEST(s21_rem_1) {
    matrix_t a = s21_create_matrix(1, 1);
    s21_remove_matrix(&a);
    ck_assert_int_eq(a.rows, FAILURE);
    ck_assert_int_eq(a.columns, FAILURE);
}
END_TEST



START_TEST(s21_eq_1) {
    matrix_t A, B, C, D;
    A = s21_create_matrix(2, 3);
    B = s21_create_matrix(2, 3);
    C = s21_create_matrix(2, 3);
    D = s21_create_matrix(1, 3);
    A.matrix[0][0] = 1.3455;
    A.matrix[0][1] = 4.2;
    A.matrix[0][2] = 3.4334;
    A.matrix[1][0] = 2.33;
    A.matrix[1][1] = 5.345;
    A.matrix[1][2] = 6.1;
    A.matrix_type = 0;
    B.matrix[0][0] = 1.3455;
    B.matrix[0][1] = 4.2;
    B.matrix[0][2] = 3.4334;
    B.matrix[1][0] = 2.33;
    B.matrix[1][1] = 5.345;
    B.matrix[1][2] = 6.1;
    B.matrix_type = 0;
    C.matrix[0][0] = 1.345;
    C.matrix[0][1] = 4.2;
    C.matrix[0][2] = 3.4334;
    C.matrix[1][0] = 2.33;
    C.matrix[1][1] = 5.345;
    C.matrix[1][2] = 6.1;
    D.matrix[0][0] = 1.35;
    D.matrix[0][1] = 4.2;
    D.matrix[0][2] = 3.44;
    D.matrix_type = 0;
    int res = s21_eq_matrix(&A, &B);
    int res1 = s21_eq_matrix(&A, &C);
    int res2 = s21_eq_matrix(&A, &D);
    ck_assert_int_eq(res, 1);
    ck_assert_int_eq(res1, 0);
    ck_assert_int_eq(res2, 0);
    s21_remove_matrix(&A);
    s21_remove_matrix(&B);
    s21_remove_matrix(&C);
    s21_remove_matrix(&D);
} END_TEST


START_TEST(s21_eq_2) {
    matrix_t a = s21_create_matrix(1, 1);
    matrix_t b = s21_create_matrix(1, 1);
    ck_assert_int_eq(s21_eq_matrix(&a, &b), SUCCESS);
    a.matrix[0][0] = 2;
    ck_assert_int_eq(s21_eq_matrix(&a, &b), FAILURE);
    b.matrix[0][0] = 2.00000001;
    ck_assert_int_eq(s21_eq_matrix(&a, &b), SUCCESS);
    b.matrix[0][0] = 2.0000002;
    ck_assert_int_eq(s21_eq_matrix(&a, &b), SUCCESS);
    matrix_t c = s21_create_matrix(2, 1);
    c.matrix[0][0] = 2;
    ck_assert_int_eq(s21_eq_matrix(&a, &c), FAILURE);
    matrix_t d = s21_create_matrix(1, 1);
    d.matrix[0][0] = 2;
    d.matrix_type = 2;
    ck_assert_int_eq(s21_eq_matrix(&a, &d), 1);
    s21_remove_matrix(&a);
    s21_remove_matrix(&b);
    s21_remove_matrix(&c);
    s21_remove_matrix(&d);
} END_TEST


START_TEST(s21_sum_1) {
    matrix_t A, B, C, D, E, F;
    A = s21_create_matrix(2, 3);
    B = s21_create_matrix(2, 3);
    C = s21_create_matrix(2, 3);
    E = s21_create_matrix(2, 2);
    A.matrix[0][0] = 1.3455;
    A.matrix[0][1] = 4.2;
    A.matrix[0][2] = 3.4334;
    A.matrix[1][0] = 2.33;
    A.matrix[1][1] = 5.345;
    A.matrix[1][2] = 6.1;
    A.matrix_type = 0;
    B.matrix[0][0] = 1.3436;
    B.matrix[0][1] = 3.23;
    B.matrix[0][2] = 3.434;
    B.matrix[1][0] = 2.3663;
    B.matrix[1][1] = 5.345;
    B.matrix[1][2] = 4.134;
    B.matrix_type = 0;
    C.matrix[0][0] = 2.6891;
    C.matrix[0][1] = 7.43;
    C.matrix[0][2] = 6.8674;
    C.matrix[1][0] = 4.6963;
    C.matrix[1][1] = 10.69;
    C.matrix[1][2] = 10.234;
    C.matrix_type = 0;
    D = s21_sum_matrix(&A, &B);
    F = s21_sum_matrix(&A, &E);
    int res = s21_eq_matrix(&C, &D);
    ck_assert_int_eq(res, 1);
    ck_assert_int_eq(D.matrix_type, 0);
    ck_assert_int_eq(F.matrix_type, 1);
    s21_remove_matrix(&A);
    s21_remove_matrix(&B);
    s21_remove_matrix(&C);
    s21_remove_matrix(&D);
    s21_remove_matrix(&E);
    s21_remove_matrix(&F);
} END_TEST


START_TEST(s21_sum_2) {
    matrix_t a = s21_create_matrix(2, 2);
    matrix_t b = s21_create_matrix(2, 2);
    a.matrix[0][0] = 1;
    a.matrix[1][1] = 1;
    matrix_t c = s21_sum_matrix(&a, &b);
    ck_assert_int_eq(c.matrix_type, IDENTITY_MATRIX);
    s21_remove_matrix(&c);

    matrix_t d = s21_create_matrix(2, 3);
    c = s21_sum_matrix(&a, &d);
    ck_assert_int_eq(c.matrix_type, INCORRECT_MATRIX);
    s21_remove_matrix(&c);

    matrix_t f = s21_create_matrix(2, 2);
    c = s21_sum_matrix(&a, &f);
    ck_assert_int_eq(c.matrix_type, 2);
    s21_remove_matrix(&a);
    s21_remove_matrix(&b);
    s21_remove_matrix(&c);
    s21_remove_matrix(&d);
    s21_remove_matrix(&f);
} END_TEST

START_TEST(s21_sub_1) {
    matrix_t A, B, C, D, E, F;
    A = s21_create_matrix(3, 3);
    B = s21_create_matrix(3, 3);
    C = s21_create_matrix(3, 3);
    E = s21_create_matrix(2, 2);
    A.matrix[0][0] = 2.0;
    A.matrix[0][1] = 0.0;
    A.matrix[0][2] = 0.0;
    A.matrix[1][0] = 0.0;
    A.matrix[1][1] = 2.0;
    A.matrix[1][2] = 0.0;
    A.matrix[2][0] = 0.0;
    A.matrix[2][1] = 0.0;
    A.matrix[2][2] = 2.0;
    A.matrix_type = 0;
    B.matrix[0][0] = 1.0;
    B.matrix[0][1] = 0.0;
    B.matrix[0][2] = 0.0;
    B.matrix[1][0] = 0.0;
    B.matrix[1][1] = 1.0;
    B.matrix[1][2] = 0.0;
    B.matrix[2][0] = 0.0;
    B.matrix[2][1] = 0.0;
    B.matrix[2][2] = 1.0;
    B.matrix_type = 0;
    C.matrix[0][0] = 1.0;
    C.matrix[0][1] = 0.0;
    C.matrix[0][2] = 0.0;
    C.matrix[1][0] = 0.0;
    C.matrix[1][1] = 1.0;
    C.matrix[1][2] = 0.0;
    C.matrix[2][0] = 0.0;
    C.matrix[2][1] = 0.0;
    C.matrix[2][2] = 1.0;
    C.matrix_type = 2;
    D = s21_sub_matrix(&A, &B);
    F = s21_sub_matrix(&A, &E);
    int res = s21_eq_matrix(&C, &D);
    ck_assert_int_eq(res, 1);
    ck_assert_int_eq(D.matrix_type, 2);
    ck_assert_int_eq(F.matrix_type, 1);
    s21_remove_matrix(&A);
    s21_remove_matrix(&B);
    s21_remove_matrix(&C);
    s21_remove_matrix(&D);
    s21_remove_matrix(&E);
    s21_remove_matrix(&F);
} END_TEST

START_TEST(s21_sub_2) {
    matrix_t A, B, C, D;
    A = s21_create_matrix(3, 3);
    B = s21_create_matrix(3, 3);
    C = s21_create_matrix(3, 3);
    A.matrix[0][0] = -1.0;
    A.matrix[0][1] = 0.0;
    A.matrix[0][2] = 0.0;
    A.matrix[1][0] = 0.0;
    A.matrix[1][1] = -1.0;
    A.matrix[1][2] = 0.0;
    A.matrix[2][0] = 0.0;
    A.matrix[2][1] = 0.0;
    A.matrix[2][2] = -1.0;
    A.matrix_type = 0;
    B.matrix[0][0] = -1.0;
    B.matrix[0][1] = 0.0;
    B.matrix[0][2] = 0.0;
    B.matrix[1][0] = 0.0;
    B.matrix[1][1] = -1.0;
    B.matrix[1][2] = 0.0;
    B.matrix[2][0] = 0.0;
    B.matrix[2][1] = 0.0;
    B.matrix[2][2] = -1.0;
    B.matrix_type = 0;
    D = s21_sub_matrix(&A, &B);
    int res = s21_eq_matrix(&C, &D);
    ck_assert_int_eq(res, 1);
    ck_assert_int_eq(D.matrix_type, 3);
    s21_remove_matrix(&A);
    s21_remove_matrix(&B);
    s21_remove_matrix(&C);
    s21_remove_matrix(&D);
} END_TEST


START_TEST(s21_sub_3) {
    matrix_t a = s21_create_matrix(2, 2);
    matrix_t b = s21_create_matrix(2, 2);
    a.matrix[0][0] = 2;
    a.matrix[1][1] = 2;
    b.matrix[0][0] = 1;
    b.matrix[1][1] = 1;
    matrix_t c;
    c = s21_sub_matrix(&a, &b);
    ck_assert_int_eq(c.matrix_type, IDENTITY_MATRIX);
    s21_remove_matrix(&c);
    matrix_t d = s21_create_matrix(2, 3);
    c = s21_sub_matrix(&a, &d);
    ck_assert_int_eq(c.matrix_type, 1);
    s21_remove_matrix(&c);
    matrix_t f = s21_create_matrix(2, 2);
    s21_remove_matrix(&f);
    f.matrix_type = INCORRECT_MATRIX;
    c = s21_sub_matrix(&a, &f);
    ck_assert_int_eq(c.matrix_type, INCORRECT_MATRIX);
    s21_remove_matrix(&a);
    s21_remove_matrix(&b);
    s21_remove_matrix(&c);
    s21_remove_matrix(&d);
} END_TEST

START_TEST(s21_mult_num_1) {
    matrix_t A, B, D;
    double x = 3.45;
    A = s21_create_matrix(2, 2);
    B = s21_create_matrix(2, 2);
    A.matrix[0][0] = 2.987654;
    A.matrix[0][1] = -5.3456;
    A.matrix[1][0] = 6.23;
    A.matrix[1][1] = -3.9292;
    A.matrix_type = 0;
    B.matrix[0][0] = 10.3074063;
    B.matrix[0][1] = -18.44232;
    B.matrix[1][0] = 21.4935;
    B.matrix[1][1] = -13.55574;
    B.matrix_type = 0;
    D = s21_mult_number(&A, x);
    int res = s21_eq_matrix(&D, &B);
    ck_assert_int_eq(res, 1);
    s21_remove_matrix(&A);
    s21_remove_matrix(&B);
    s21_remove_matrix(&D);
} END_TEST

START_TEST(s21_mult_num_2) {
    matrix_t a = s21_create_matrix(2, 2);
    matrix_t c;
    a.matrix[0][0] = 1;
    a.matrix[1][1] = 1;
    c = s21_mult_number(&a, 1);
    ck_assert_int_eq(c.matrix_type, IDENTITY_MATRIX);
    s21_remove_matrix(&c);
    c = s21_mult_number(&a, 2);
    int check = 0;
    for (int i = 0; i < c.rows; i++) {
        for (int j = 0; j < c.columns; j++) {
            if (c.matrix[i][j] != (a.matrix[i][j] * 2)) {
                check = 1;
                break;
            }
            if (check) {
                break;
            }
        }
    }
    s21_remove_matrix(&c);
    ck_assert_int_eq(check, FAILURE);
    s21_remove_matrix(&a);
    a.matrix_type = INCORRECT_MATRIX;
    c = s21_mult_number(&a, 2);
    ck_assert_int_eq(c.matrix_type, INCORRECT_MATRIX);
    s21_remove_matrix(&c);
} END_TEST

START_TEST(s21_mult_m_1) {
    matrix_t A, B, C, D, E, F;
    A = s21_create_matrix(2, 2);
    B = s21_create_matrix(2, 4);
    C = s21_create_matrix(2, 4);
    E = s21_create_matrix(4, 2);
    A.matrix[0][0] = 2.987654;
    A.matrix[0][1] = -5.3456;
    A.matrix[1][0] = 6.23;
    A.matrix[1][1] = -3.9292;
    A.matrix_type = 0;
    B.matrix[0][0] = 4.838;
    B.matrix[0][1] = 8.0987;
    B.matrix[0][2] = 8.09;
    B.matrix[0][3] = -987;
    B.matrix[1][0] = -8;
    B.matrix[1][1] = -0.921;
    B.matrix[1][2] = -0.91;
    B.matrix[1][3] = 21;
    B.matrix_type = 0;
    C.matrix[0][0] = 57.219070052;
    C.matrix[0][1] = 29.1194110498;
    C.matrix[0][2] = 29.03461686;
    C.matrix[0][3] = -3061.072098;
    C.matrix[1][0] = 61.57434;
    C.matrix[1][1] = 54.0736942;
    C.matrix[1][2] = 53.976272;
    C.matrix[1][3] = -6231.5232;
    C.matrix_type = 0;
    D = s21_mult_matrix(&A, &B);
    F = s21_mult_matrix(&A, &E);
    int res = s21_eq_matrix(&C, &D);
    ck_assert_int_eq(res, 1);
    ck_assert_int_eq(C.matrix_type, 0);
    ck_assert_int_eq(F.matrix_type, 1);
    s21_remove_matrix(&A);
    s21_remove_matrix(&B);
    s21_remove_matrix(&C);
    s21_remove_matrix(&D);
    s21_remove_matrix(&E);
    s21_remove_matrix(&F);
} END_TEST


START_TEST(s21_mult_m_2) {
    matrix_t a = s21_create_matrix(3, 2);
    matrix_t b = s21_create_matrix(2, 3);
    a.matrix[0][0] = 1;
    a.matrix[0][1] = 4;
    a.matrix[1][0] = 2;
    a.matrix[1][1] = 5;
    a.matrix[2][0] = 3;
    a.matrix[2][1] = 6;
    b.matrix[0][0] = 1;
    b.matrix[0][1] = -1;
    b.matrix[0][2] = 1;
    b.matrix[1][0] = 2;
    b.matrix[1][1] = 3;
    b.matrix[1][2] = 4;
    matrix_t res = s21_create_matrix(3, 3);
    res.matrix[0][0] = 9;
    res.matrix[0][1] = 11;
    res.matrix[0][2] = 17;
    res.matrix[1][0] = 12;
    res.matrix[1][1] = 13;
    res.matrix[1][2] = 22;
    res.matrix[2][0] = 15;
    res.matrix[2][1] = 15;
    res.matrix[2][2] = 27;
    matrix_t c = s21_mult_matrix(&a, &b);
    int check = 0;
    for (int i = 0; i < c.rows; i++) {
        for (int j = 0; j < c.columns; j++) {
            if (c.matrix[i][j] != res.matrix[i][j]) {
                check = 1;
                break;
            }
            if (check) {
                break;
            }
        }
    }
    s21_remove_matrix(&c);
    ck_assert_int_eq(check, FAILURE);
    s21_remove_matrix(&a);
    a.matrix_type = INCORRECT_MATRIX;
    c = s21_mult_matrix(&a, &b);
    ck_assert_int_eq(c.matrix_type, INCORRECT_MATRIX);
    s21_remove_matrix(&c);
    s21_remove_matrix(&b);
    s21_remove_matrix(&res);
    matrix_t m1 = s21_create_matrix(3, 1);
    s21_remove_matrix(&m1);
    m1.matrix_type = INCORRECT_MATRIX;
    matrix_t m2 = s21_create_matrix(2, 4);
    matrix_t m3 = s21_mult_matrix(&m1, &m2);
    ck_assert_int_eq(m3.matrix_type, INCORRECT_MATRIX);
    s21_remove_matrix(&m3);
    s21_remove_matrix(&m2);
} END_TEST


START_TEST(s21_trans_1) {
    matrix_t A, B, C;
    A = s21_create_matrix(3, 3);
    C = s21_create_matrix(3, 3);
    A.matrix[0][0] = 2.0;
    A.matrix[0][1] = 0.0;
    A.matrix[0][2] = 0.0;
    A.matrix[1][0] = 0.0;
    A.matrix[1][1] = 2.0;
    A.matrix[1][2] = 0.0;
    A.matrix[2][0] = 0.0;
    A.matrix[2][1] = 0.0;
    A.matrix[2][2] = 2.0;
    A.matrix_type = 0;
    B = s21_transpose(&A);
    C.matrix[0][0] = 2.0;
    C.matrix[0][1] = 0.0;
    C.matrix[0][2] = 0.0;
    C.matrix[1][0] = 0.0;
    C.matrix[1][1] = 2.0;
    C.matrix[1][2] = 0.0;
    C.matrix[2][0] = 0.0;
    C.matrix[2][1] = 0.0;
    C.matrix[2][2] = 2.0;
    C.matrix_type = 0;
    int res = s21_eq_matrix(&B, &C);
    ck_assert_int_eq(res, 1);
    s21_remove_matrix(&A);
    s21_remove_matrix(&B);
    s21_remove_matrix(&C);
} END_TEST

START_TEST(s21_trans_2) {
    matrix_t a = s21_create_matrix(2, 1);
    matrix_t b = s21_create_matrix(1, 2);
    a.matrix[0][0] = 1;
    a.matrix[1][0] = 2;
    a.matrix_type = 0;
    b.matrix[0][0] = 1;
    b.matrix[0][1] = 2;
    b.matrix_type = 0;
    matrix_t c = s21_transpose(&a);
    ck_assert_int_eq(s21_eq_matrix(&b, &c), SUCCESS);
    s21_remove_matrix(&c);
    s21_remove_matrix(&a);
    a.matrix_type = INCORRECT_MATRIX;
    c = s21_transpose(&a);
    ck_assert_int_eq(s21_eq_matrix(&b, &c), FAILURE);
    s21_remove_matrix(&c);
    s21_remove_matrix(&b);
} END_TEST

START_TEST(s21_calc_compl_1) {
    matrix_t A, B, C, D, E;
    A = s21_create_matrix(3, 3);
    C = s21_create_matrix(3, 3);
    D = s21_create_matrix(2, 3);
    A.matrix[0][0] = 1.0;
    A.matrix[0][1] = 2.0;
    A.matrix[0][2] = 3.0;
    A.matrix[1][0] = 0.0;
    A.matrix[1][1] = 4.0;
    A.matrix[1][2] = 2.0;
    A.matrix[2][0] = 5.0;
    A.matrix[2][1] = 2.0;
    A.matrix[2][2] = 1.0;
    A.matrix_type = 0;
    B = s21_calc_complements(&A);
    C.matrix[0][0] = 0.0;
    C.matrix[0][1] = 10.0;
    C.matrix[0][2] = -20.0;
    C.matrix[1][0] = 4.0;
    C.matrix[1][1] = -14.0;
    C.matrix[1][2] = 8.0;
    C.matrix[2][0] = -8.0;
    C.matrix[2][1] = -2.0;
    C.matrix[2][2] = 4.0;
    C.matrix_type = 0;
    int res = s21_eq_matrix(&B, &C);
    ck_assert_int_eq(res, 1);
    ck_assert_int_eq(B.matrix_type, 0);
    E = s21_calc_complements(&D);
    ck_assert_int_eq(E.matrix_type, 1);
    s21_remove_matrix(&A);
    s21_remove_matrix(&B);
    s21_remove_matrix(&C);
    s21_remove_matrix(&D);
    s21_remove_matrix(&E);
} END_TEST

START_TEST(s21_calc_compl_2) {
    matrix_t a = s21_create_matrix(3, 3);
    matrix_t res = s21_create_matrix(3, 3);
    matrix_t m;
    res.matrix[0][0] = 0;
    res.matrix[0][1] = 10;
    res.matrix[0][2] = -20;
    res.matrix[1][0] = 4;
    res.matrix[1][1] = -14;
    res.matrix[1][2] = 8;
    res.matrix[2][0] = -8;
    res.matrix[2][1] = -2;
    res.matrix[2][2] = 4;
    res.matrix_type = 0;
    a.matrix[0][0] = 1;
    a.matrix[0][1] = 2;
    a.matrix[0][2] = 3;
    a.matrix[1][0] = 0;
    a.matrix[1][1] = 4;
    a.matrix[1][2] = 2;
    a.matrix[2][0] = 5;
    a.matrix[2][1] = 2;
    a.matrix[2][2] = 1;
    a.matrix_type = 0;
    m = s21_calc_complements(&a);
    ck_assert_int_eq(s21_eq_matrix(&res, &m), SUCCESS);
    s21_remove_matrix(&m);
    s21_remove_matrix(&a);
    a.matrix_type = INCORRECT_MATRIX;
    m = s21_calc_complements(&a);
    ck_assert_int_eq(s21_eq_matrix(&res, &m), FAILURE);
    s21_remove_matrix(&m);
    s21_remove_matrix(&res);
} END_TEST

START_TEST(s21_determ_1) {
    matrix_t A, B, C, D;
    A = s21_create_matrix(2, 2);
    B = s21_create_matrix(2, 4);
    C = s21_create_matrix(1, 1);
    D = s21_create_matrix(3, 3);
    A.matrix[0][0] = 2.987654;
    A.matrix[0][1] = -5.3456;
    A.matrix[1][0] = 6.23;
    A.matrix[1][1] = -3.9292;
    A.matrix_type = 0;
    B.matrix[0][0] = 4.838;
    B.matrix[0][1] = 8.0987;
    B.matrix[0][2] = 8.09;
    B.matrix[0][3] = -987;
    B.matrix[1][0] = -8;
    B.matrix[1][1] = -0.921;
    B.matrix[1][2] = -0.91;
    B.matrix[1][3] = 21;
    B.matrix_type = 0;
    C.matrix[0][0] = 2.32443;
    D.matrix[0][0] = 2.8;
    D.matrix[0][1] = 1.3;
    D.matrix[0][2] = 7.01;
    D.matrix[1][0] = -1.03;
    D.matrix[1][1] = -2.3;
    D.matrix[1][2] = 3.01;
    D.matrix[2][0] = 0;
    D.matrix[2][1] = -3;
    D.matrix[2][2] = 2;
    double res = s21_determinant(&A);
    ck_assert_float_eq(res, 21.563997);
    double res1 = s21_determinant(&C);
    ck_assert_float_eq(res1, 2.32443);
    s21_determinant(&B);
    double res2 = s21_determinant(&D);
    ck_assert_float_eq(res2, 36.7429);
    s21_remove_matrix(&A);
    s21_remove_matrix(&B);
    s21_remove_matrix(&C);
    s21_remove_matrix(&D);
} END_TEST

START_TEST(s21_determ_2) {
    matrix_t d = s21_create_matrix(4, 4);
    d.matrix[0][0] = 2;
    d.matrix[0][1] = 5;
    d.matrix[0][2] = 8;
    d.matrix[0][3] = 10;
    d.matrix[1][0] = 2;
    d.matrix[1][1] = 4;
    d.matrix[1][2] = 7;
    d.matrix[1][3] = 6;
    d.matrix[2][0] = 5;
    d.matrix[2][1] = 1;
    d.matrix[2][2] = 2;
    d.matrix[2][3] = 3;
    d.matrix[3][0] = 5;
    d.matrix[3][1] = 8;
    d.matrix[3][2] = 9;
    d.matrix[3][3] = 0;
    double det = s21_determinant(&d);
    ck_assert_double_eq(det, -403.0);
    matrix_t a = s21_create_matrix(3, 3);
    a.matrix[0][0] = 1;
    a.matrix[0][1] = 2;
    a.matrix[0][2] = 3;
    a.matrix[1][0] = 4;
    a.matrix[1][1] = 5;
    a.matrix[1][2] = 6;
    a.matrix[2][0] = 7;
    a.matrix[2][1] = 8;
    a.matrix[2][2] = 9;
    det = s21_determinant(&a);
    ck_assert_double_eq(det, 0.0);
    s21_remove_matrix(&a);
    a.matrix_type = INCORRECT_MATRIX;
    det = s21_determinant(&a);
    ck_assert_msg(isnan(det), "determinant is nan");
    s21_remove_matrix(&d);
} END_TEST

START_TEST(s21_inverse_1) {
    matrix_t A, C;
    A = s21_create_matrix(3, 3);
    C = s21_create_matrix(3, 3);
    C.matrix[0][0] = 1.0;
    C.matrix[0][1] = -1.0;
    C.matrix[0][2] = 1.0;
    C.matrix[1][0] = -38.0;
    C.matrix[1][1] = 41.0;
    C.matrix[1][2] = -34.0;
    C.matrix[2][0] = 27.0;
    C.matrix[2][1] = -29.0;
    C.matrix[2][2] = 24.0;
    C.matrix_type = 0;
    A.matrix[0][0] = 2.0;
    A.matrix[0][1] = 5.0;
    A.matrix[0][2] = 7.0;
    A.matrix[1][0] = 6.0;
    A.matrix[1][1] = 3.0;
    A.matrix[1][2] = 4.0;
    A.matrix[2][0] = 5.0;
    A.matrix[2][1] = -2.0;
    A.matrix[2][2] = -3.0;
    A.matrix_type = 0;
    matrix_t B = s21_inverse_matrix(&A);
    int res = s21_eq_matrix(&B, &C);
    ck_assert_int_eq(res, 1);
    s21_remove_matrix(&A);
    s21_remove_matrix(&B);
    s21_remove_matrix(&C);
} END_TEST

START_TEST(s21_inverse_2) {
    matrix_t A, B, C, D, E, F, G;
    A = s21_create_matrix(3, 3);
    C = s21_create_matrix(3, 3);
    D = s21_create_matrix(2, 3);
    F = s21_create_matrix(2, 2);
    C.matrix[0][0] = 44300.0/367429.0;
    C.matrix[0][1] = -236300.0/367429.0;
    C.matrix[0][2] = 200360.0/367429.0;
    C.matrix[1][0] = 20600.0/367429.0;
    C.matrix[1][1] = 56000.0/367429.0;
    C.matrix[1][2] = -156483.0/367429.0;
    C.matrix[2][0] = 30900.0/367429.0;
    C.matrix[2][1] = 84000.0/367429.0;
    C.matrix[2][2] = -51010.0/367429.0;
    C.matrix_type = 0;
    A.matrix[0][0] = 2.8;
    A.matrix[0][1] = 1.3;
    A.matrix[0][2] = 7.01;
    A.matrix[1][0] = -1.03;
    A.matrix[1][1] = -2.3;
    A.matrix[1][2] = 3.01;
    A.matrix[2][0] = 0;
    A.matrix[2][1] = -3;
    A.matrix[2][2] = 2;
    A.matrix_type = 0;
    B = s21_inverse_matrix(&A);
    int res = s21_eq_matrix(&B, &C);
    ck_assert_int_eq(res, 1);
    E = s21_inverse_matrix(&D);
    ck_assert_int_eq(E.matrix_type, 1);
    G = s21_inverse_matrix(&F);
    ck_assert_int_eq(G.matrix_type, 1);
    s21_remove_matrix(&A);
    s21_remove_matrix(&B);
    s21_remove_matrix(&C);
    s21_remove_matrix(&D);
    s21_remove_matrix(&F);
} END_TEST


START_TEST(s21_inverse_3) {
    matrix_t res = s21_create_matrix(3, 3);
    res.matrix[0][0] = 1;
    res.matrix[0][1] = -2;
    res.matrix[0][2] = 1;
    res.matrix[1][0] = -2;
    res.matrix[1][1] = 2;
    res.matrix[1][2] = -0.5;
    res.matrix[2][0] = 1;
    res.matrix[2][1] = -0.5;
    res.matrix[2][2] = 0;
    res.matrix_type = 0;

    matrix_t a = s21_create_matrix(3, 3);
    a.matrix[0][0] = 1;
    a.matrix[0][1] = 2;
    a.matrix[0][2] = 4;
    a.matrix[1][0] = 2;
    a.matrix[1][1] = 4;
    a.matrix[1][2] = 6;
    a.matrix[2][0] = 4;
    a.matrix[2][1] = 6;
    a.matrix[2][2] = 8;
    a.matrix_type = 0;

    matrix_t a_obr = s21_inverse_matrix(&a);

    ck_assert_int_eq(s21_eq_matrix(&res, &a_obr), SUCCESS);
    ck_assert_int_eq(a_obr.matrix_type, CORRECT_MATRIX);
    s21_remove_matrix(&a_obr);

    a.matrix_type = INCORRECT_MATRIX;
    a_obr = s21_inverse_matrix(&a);
    ck_assert_int_eq(a_obr.matrix_type, INCORRECT_MATRIX);

    s21_remove_matrix(&a);
    s21_remove_matrix(&res);
} END_TEST

int main(void) {
  Suite *s1 = suite_create("s21_matrix");
  TCase *tc1_1 = tcase_create("s21_matrix");
  SRunner *sr = srunner_create(s1);
  int nf;
  suite_add_tcase(s1, tc1_1);
  tcase_add_test(tc1_1, s21_create_1);
  tcase_add_test(tc1_1, s21_create_2);
  tcase_add_test(tc1_1, s21_create_3);
  tcase_add_test(tc1_1, s21_rem_1);
  tcase_add_test(tc1_1, s21_eq_1);
  tcase_add_test(tc1_1, s21_eq_2);
  tcase_add_test(tc1_1, s21_sum_1);
  tcase_add_test(tc1_1, s21_sum_2);
  tcase_add_test(tc1_1, s21_sub_1);
  tcase_add_test(tc1_1, s21_sub_2);
  tcase_add_test(tc1_1, s21_sub_3);
  tcase_add_test(tc1_1, s21_mult_num_1);
  tcase_add_test(tc1_1, s21_mult_num_2);
  tcase_add_test(tc1_1, s21_mult_m_1);
  tcase_add_test(tc1_1, s21_mult_m_2);
  tcase_add_test(tc1_1, s21_trans_1);
  tcase_add_test(tc1_1, s21_trans_2);
  tcase_add_test(tc1_1, s21_calc_compl_1);
  tcase_add_test(tc1_1, s21_calc_compl_2);
  tcase_add_test(tc1_1, s21_determ_1);
  tcase_add_test(tc1_1, s21_determ_2);
  tcase_add_test(tc1_1, s21_inverse_1);
  tcase_add_test(tc1_1, s21_inverse_2);
  tcase_add_test(tc1_1, s21_inverse_3);
  srunner_run_all(sr, CK_ENV);
  nf = srunner_ntests_failed(sr);
  srunner_free(sr);
  return nf == 0 ? 0 : 1;
}
