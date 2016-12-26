/*
	2016/07/14
	カタツムリの移動シミュレーション
	藤井涼平
*/

#include <math.h>
#include <stdio.h>

//実行ステップ数
#define STEPS 10000
//時間刻み幅
#define DELTA 0.01

//質点数
#define POINTS 50
//質点質量
#define M 0.0075

//ダンパー粘性係数
#define Q 5e-3
//バネ定数の比例係数
#define KAPPA 1.0
//平均自然長
#define L 1

//波数
#define BETA 0.16
//位相速度
#define OMEGA 0.8
//振幅
#define ALPHA 0.24

//粘液降伏値
#define FU 3e-2
//粘液回復値
#define FL 1e-3
//粘液（弾性固体）に働く力の比例定数
#define GAMMA 1e-3
//粘液（粘性流体）に働く力の比例定数
#define ETA 1e-6
////粘液（弾性固体）に働く力の比例定数
//#define GAMMA 0
////粘液（粘性流体）に働く力の比例定数
//#define ETA 0


/**
初期化
*/
void initialize(void);
/**
位置を文字列として表示
@param time 時刻
*/
void toString(double time);
/**
4次Runge-Kutta法を用いて位置を更新
@param currentTime 現在時刻
*/
void update(double currentTime);
/**
ベクトルを加算
@param vec1, vec2 加算するベクトル
@param result 結果を格納するベクトル
@return 結果を格納するベクトルへのポインタ
*/
double* add(double vec1[], double vec2[], double result[]);
/**
ベクトルを加算
@param vec1, vec2 減算するベクトル（vec1 - vec2）
@param result 結果を格納するベクトル
@return 結果を格納するベクトルへのポインタ
*/
double* substract(double vec1[], double vec2[], double result[]);
/**
ベクトルとスカラーを積算
@param vec1 ベクトル
@param value スカラー
@param result 結果を格納するベクトル
@return 結果を格納するベクトルへのポインタ
*/
double* multiply(double vec1[], double value, double result[]);
/**
バネの自然長
@param n バネの番号
@param time 時刻
@return 自然長
*/
double l(int n, double time);
/**
バネ定数
@param n バネの番号
@param time 時刻
@return 自然長
*/
double k(int n, double time);
/**
粘液なしの加速度
@param n バネの番号
@param time 時刻
@param x[] 位置の配列
@param v[] 速度の配列
@return 粘液なしの加速度
*/
double a_nomucus(int n, double time, double x[], double v[]);
/**
粘液ありの加速度
@param n バネの番号
@param time 時刻
@param x[] 位置の配列
@param v[] 速度の配列
@return 粘液ありの加速度
*/
double a(int n, double time, double x[], double v[]);

//位置配列
double x_g[POINTS];

//速度配列
double v_g[POINTS];

//最後に弾性固体になった位置の配列
double xBar[POINTS];

//弾性固体(1)／粘性流体(0)
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
	//Runge-Kuttaで用いるkベクトル
	double kx1[POINTS];
	double kx2[POINTS];
	double kx3[POINTS];
	double kx4[POINTS];
	double kv1[POINTS];
	double kv2[POINTS];
	double kv3[POINTS];
	double kv4[POINTS];

	//ベクトルの計算結果を一時的に格納する配列
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