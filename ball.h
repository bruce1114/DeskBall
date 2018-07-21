#ifndef BALL
#define BALL
#include<graphics.h>
typedef struct ball {
	int x, y, xl, yl;//当前，上一次坐标
	int r;//半径
	bool move;//是否运动
	bool movech;//。。。
	bool isexist;//是否存在，掉洞不存在
	double v;//速度
	double angle;//方向角
	double morex;//x方向位移累积
	double morey;//同上
	double mass;//质量
	ball();
}ball;

typedef struct hole {
	int x, y;
	int r;
	bool repaint;
	hole();
}hole;

typedef struct Pair {
	int x[2], y[2];
}Pair;

typedef struct zhangai {
	int pairnum;
	Pair p[10];
}zhangai;

typedef struct cirpoint {
	int x,y;
}cirpoint;

extern double Time;
extern double PI;
extern double ug;
extern int cnt;
extern int kuan, gao;
extern int extra;
extern int inballnum;
extern int sleeptime;
extern MOUSEMSG m;
extern int mlx, mly;
extern double maxlen;
extern double rate;
extern int mainball;//const
extern int holenum;
extern int ballnum;
extern bool self;
extern ball test[];
extern hole h[];
extern int bkcolor, mainbcolor, otherbcolor, holecolor, linecolor, zhangaicolor;
extern bool nocancel;
extern int occupy[2000][1000];
extern cirpoint cirp[];
extern zhangai zhang[];
extern int zhangnum;
extern void(*builder)();
extern void(*drawer)();
extern int bigball;
extern int solidball;
extern int ch;

extern double just;

void move(ball&,double);
void cal(ball*, int);
void rebound(ball&, ball&,double);
void calvecangle(double&, double, double);
bool inhole(ball*, int, hole*, int);
void repaint(hole*, int);
void mouseevent(ball&);
void inibuilder(int*, int*,int);
void builder0();
void drawer0();
void iszhangai(zhangai*,int,ball*,int);
bool rebound1(cirpoint&, ball&);
void rebound2(int, ball&);
bool judgepair(Pair&, ball&);
bool judgepair1(Pair&, ball&);
void choice();
void drawer1();
void builder1();
void drawer2();
void builder2();
void drawer3();
void builder3();
void choice5();
void drawern();
void buildern();
#endif // !BALL
