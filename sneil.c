/*
	2016/07/14
	�J�^�c�����̈ړ��V�~�����[�V����
	�������
*/

#include <math.h>
#include <stdio.h>

//���s�X�e�b�v��
#define STEPS 10000
//���ԍ��ݕ�
#define DELTA 0.01

//���_��
#define POINTS 50
//���_����
#define M 0.0075

//�_���p�[�S���W��
#define Q 5e-3
//�o�l�萔�̔��W��
#define KAPPA 1.0
//���ώ��R��
#define L 1

//�g��
#define BETA 0.16
//�ʑ����x
#define OMEGA 0.8
//�U��
#define ALPHA 0.24

//�S�t�~���l
#define FU 3e-2
//�S�t�񕜒l
#define FL 1e-3
//�S�t�i�e���ő́j�ɓ����͂̔��萔
#define GAMMA 1e-3
//�S�t�i�S�����́j�ɓ����͂̔��萔
#define ETA 1e-6
////�S�t�i�e���ő́j�ɓ����͂̔��萔
//#define GAMMA 0
////�S�t�i�S�����́j�ɓ����͂̔��萔
//#define ETA 0


/**
������
*/
void initialize(void);
/**
�ʒu�𕶎���Ƃ��ĕ\��
@param time ����
*/
void toString(double time);
/**
4��Runge-Kutta�@��p���Ĉʒu���X�V
@param currentTime ���ݎ���
*/
void update(double currentTime);
/**
�x�N�g�������Z
@param vec1, vec2 ���Z����x�N�g��
@param result ���ʂ��i�[����x�N�g��
@return ���ʂ��i�[����x�N�g���ւ̃|�C���^
*/
double* add(double vec1[], double vec2[], double result[]);
/**
�x�N�g�������Z
@param vec1, vec2 ���Z����x�N�g���ivec1 - vec2�j
@param result ���ʂ��i�[����x�N�g��
@return ���ʂ��i�[����x�N�g���ւ̃|�C���^
*/
double* substract(double vec1[], double vec2[], double result[]);
/**
�x�N�g���ƃX�J���[��ώZ
@param vec1 �x�N�g��
@param value �X�J���[
@param result ���ʂ��i�[����x�N�g��
@return ���ʂ��i�[����x�N�g���ւ̃|�C���^
*/
double* multiply(double vec1[], double value, double result[]);
/**
�o�l�̎��R��
@param n �o�l�̔ԍ�
@param time ����
@return ���R��
*/
double l(int n, double time);
/**
�o�l�萔
@param n �o�l�̔ԍ�
@param time ����
@return ���R��
*/
double k(int n, double time);
/**
�S�t�Ȃ��̉����x
@param n �o�l�̔ԍ�
@param time ����
@param x[] �ʒu�̔z��
@param v[] ���x�̔z��
@return �S�t�Ȃ��̉����x
*/
double a_nomucus(int n, double time, double x[], double v[]);
/**
�S�t����̉����x
@param n �o�l�̔ԍ�
@param time ����
@param x[] �ʒu�̔z��
@param v[] ���x�̔z��
@return �S�t����̉����x
*/
double a(int n, double time, double x[], double v[]);

//�ʒu�z��
double x_g[POINTS];

//���x�z��
double v_g[POINTS];

//�Ō�ɒe���ő̂ɂȂ����ʒu�̔z��
double xBar[POINTS];

//�e���ő�(1)�^�S������(0)
int sigma[POINTS];

int main(void) {
	int i;
	initialize();

	toString(0);
	for (i = 0; i < STEPS; i++) {
		update(DELTA * i);
		if (i % (STEPS / 500) == 0)
			toString(DELTA * i);
	}
	return 0;
}

void initialize(void) {
	int i;
	double temp = 0;
	for (i = 0; i < POINTS; i++) {
		v_g[i] = 0;
		temp += l(i, 0);
		x_g[i] = temp;
		xBar[i] = temp;
		sigma[i] = 1;
	}
}

void toString(double time) {
	int i;
	for (i = 0; i < POINTS; i++) {
		printf("%f\t%f\n", time, x_g[i]);
	}
}

void update(double currentTime) {
	//Runge-Kutta�ŗp����k�x�N�g��
	double kx1[POINTS];
	double kx2[POINTS];
	double kx3[POINTS];
	double kx4[POINTS];
	double kv1[POINTS];
	double kv2[POINTS];
	double kv3[POINTS];
	double kv4[POINTS];

	//�x�N�g���̌v�Z���ʂ��ꎞ�I�Ɋi�[����z��
	double temp1[POINTS];
	double temp2[POINTS];

	int i;

	for (i = 0; i < POINTS; i++) {
		kv1[i] = DELTA * a(i, currentTime, x_g, v_g);
		kx1[i] = DELTA * v_g[i];
	}

	for (i = 0; i < POINTS; i++) {
		kv2[i] = DELTA *
			a(i, currentTime + DELTA / 2,
				add(x_g, multiply(kx1, 0.5, temp1), temp2),
				add(v_g, multiply(kv1, 0.5, temp1), temp2));

		kx2[i] = DELTA *(v_g[i] + 0.5*kv1[i]);
	}

	for (i = 0; i < POINTS; i++) {
		kv3[i] = DELTA *
			a(i, currentTime + DELTA / 2,
				add(x_g, multiply(kx2, 0.5, temp1), temp2),
				add(v_g, multiply(kv2, 0.5, temp1), temp2));

		kx3[i] = DELTA *(v_g[i] + 0.5*kv2[i]);
	}

	for (i = 0; i < POINTS; i++) {
		kv4[i] = DELTA *
			a(i, currentTime + DELTA,
				add(x_g, kx3, temp1),
				add(v_g, kv3, temp1));

		kx4[i] = DELTA *(v_g[i] + kv3[i]);
	}

	for (i = 0; i < POINTS; i++) {
		v_g[i] += (kv1[i] + 2 * kv2[i] + 2 * kv3[i] + kv4[i]) / 6;
		x_g[i] += (kx1[i] + 2 * kx2[i] + 2 * kx3[i] + kx4[i]) / 6;
	}
}

double* add(double vec1[], double vec2[], double result[]) {
	int i;

	for (i = 0; i < POINTS; i++) {
		result[i] = vec1[i] + vec2[i];
	}

	return result;
}

double* substract(double vec1[], double vec2[], double result[]) {
	int i;

	for (i = 0; i < POINTS; i++) {
		result[i] = vec1[i] - vec2[i];
	}

	return result;
}

double* multiply(double vec1[], double value, double result[]) {
	int i;

	for (i = 0; i < POINTS; i++) {
		result[i] = vec1[i] * value;
	}

	return result;
}

double l(int n, double time) {
	return L*(1 + ALPHA*sin(OMEGA*time - BETA*n));
}

double k(int n, double time) {
	double length = l(n, time);

	return KAPPA / length;
}

double a_nomucus(int n, double time, double x[], double v[]) {
	double comp1;
	double comp2;
	double comp3;
	double comp4;

	if (n < 0 || n >= POINTS)
		return NAN;

	if (n == 0) {
		comp1 = k(n, time)*(x[n + 1] - x[n] - l(n, time)) / M;
		comp2 = 0;
		comp3 = Q*(v[n + 1] - v[n]);
		comp4 = 0;

		return comp1 + comp2 + comp3 + comp4;
	}

	if (n == POINTS - 1) {
		comp1 = 0;
		comp2 = -k(n - 1, time)*(x[n] - x[n - 1] - l(n - 1, time)) / M;
		comp3 = 0;
		comp4 = -Q*(v[n] - v[n - 1]);

		return comp1 + comp2 + comp3 + comp4;
	}

	comp1 = k(n, time)*(x[n + 1] - x[n] - l(n, time)) / M;
	comp2 = -k(n - 1, time)*(x[n] - x[n - 1] - l(n - 1, time)) / M;
	comp3 = Q*(v[n + 1] - v[n]);
	comp4 = -Q*(v[n] - v[n - 1]);

	return comp1 + comp2 + comp3 + comp4;
}

double a(int n, double time, double x[], double v[]) {
	double comp1234 = a_nomucus(n, time, x, v);
	double comp5;
	double comp6;
	double force = M*comp1234;

	if (fabs(force) >= FU && sigma[n] == 1) {
		sigma[n] = 0;
	}

	else if (fabs(force) <= FL && sigma[n] == 0) {
		sigma[n] = 1;
		xBar[n] = x_g[n];
	}

	comp5 = -(1 - sigma[n])*ETA*v[n]/M;
	comp6 = -sigma[n]*GAMMA*(x[n] - xBar[n])/M;

	return comp1234 + comp5 + comp6;
}