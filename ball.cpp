#include<cmath>
#include<graphics.h>
#include<iostream>
#include<cmath>
#include<conio.h>
#include"ball.h"
using namespace std;
double Time = 1;//
double PI = acos(-1);
double ug = 0.07;//
double Exp = 1e-7;
int sleeptime = 19;//19
int kuan = 1000, gao = 500;
MOUSEMSG m;
int mlx=-10, mly=-10;
double maxlen = 200.0;
double rate = 10.0;
int mainball = 0;
int bkcolor = DARKGRAY, mainbcolor = WHITE, otherbcolor = YELLOW, holecolor = BLACK, linecolor = YELLOW, zhangaicolor = WHITE;
bool nocancel = 1;
ball test[500];
hole h[20];
int occupy[2000][1000];
void(*builder)();
void(*drawer)();
cirpoint cirp[25];
zhangai zhang[20];
int holenum;
int ballnum;
int zhangnum;
int br = 11;
int brbr = br*br;
int houdu = 50;
int extra = 0;
int inballnum = 0;
int bigball = 0;
int solidball = 0;
int ch;
bool self = 0;

double just = 0;

void move(ball &test,double sTime)
{
	if (test.move&&test.isexist)
	{
		if (test.movech)
		{
			double ugTime = ug*sTime;
			double tempx;
			if (test.v <= ugTime) tempx = test.v*test.v / (ug + ug);
			else tempx = test.v*sTime - 0.5*ugTime*sTime;//
			test.morex += tempx*cos(test.angle);
			double temp = 0.0;
			if (test.morex >= 1.0)
			{
				temp = floor(test.morex);
				test.x += temp;
				test.morex -= temp;
			}
			else if (test.morex <= -1.0)
			{
				temp = ceil(test.morex);
				test.x += temp;
				test.morex -= temp;
			}
			test.morey += tempx*sin(test.angle);
			if (test.morey >= 1.0)
			{
				temp = floor(test.morey);
				test.y += temp;
				test.morey -= temp;
			}
			else if (test.morey <= -1.0)
			{
				temp = ceil(test.morey);
				test.y += temp;
				test.morey -= temp;
			}
			test.v -= ugTime;
			if (test.v <= 0)
				test.v = 0, test.move = 0;
			if (test.x + test.r > kuan)//
			{
				test.x = kuan + kuan - test.x - test.r - test.r;
				test.angle = PI - test.angle;
			}
			else if (test.x - test.r < 0)
			{
				test.x = test.r + test.r - test.x;
				test.angle = PI - test.angle;
			}
			if (test.y + test.r > gao)
			{
				test.y = gao + gao - test.y - test.r - test.r;
				test.angle = -test.angle;
			}
			else if (test.y - test.r < 0)
			{
				test.y = test.r + test.r - test.y;
				test.angle = -test.angle;
			}
		}
		else if(fabs(Time-sTime)<=Exp)
		{
			test.movech = 1;
		}
	}
}

ball::ball()
{
	move = 1;
	movech = 1;
	isexist = 1;
	r = 11;
	v = 0;
	angle = PI/3;
	morex = morey = 0;
	mass = 5;
	x = 200;
	y = 20;
}

void calvecangle(double &angle, double vecx, double vecy)
{
	angle = atan(vecy / vecx);
	if (fabs(vecx - 0.0) <= Exp&&fabs(vecy - 0.0) <= Exp) angle = 0.0;
	else if (fabs(angle - 0.0) <= Exp) angle = (vecx > 0.0 ? 0.0 : PI);
	else if (angle > 0.0) angle = (vecy > 0.0 ? angle : angle - PI);
	else if (angle < 0.0) angle = (vecy > 0.0 ? angle + PI : angle);
}

void rebound(ball &a, ball &b,double more)
{
	//计算球心到球心的向量
	double vecy = b.y - a.y;
	double vecx = b.x - a.x;
	double alfa = 0;//坐标变换的基角
	double part2 = a.mass + b.mass;
	calvecangle(alfa, vecx, vecy);//根据向量计算基角
	if (fabs(a.v - 0.0) <= Exp&&fabs(b.v - 0.0) <= Exp) goto skip;//若两球静止但意外被判相碰，就单纯的将它们推开即可
	double anglea = a.angle - alfa;//变换坐标系
	double angleb = b.angle - alfa;
	double vax = a.v*cos(anglea), vay = a.v*sin(anglea), vbx = b.v*cos(angleb), vby = b.v*sin(angleb);
	if (vax - vbx <= 0.0) goto skip;//两球被判相碰但并不满足相碰条件。。。
	double part1 = vax - vbx;
	double axv = (a.mass - b.mass) / part2*part1 + vbx;
	double bxv = (a.mass + a.mass) / part2*part1 + vbx;
	double angleaN = 0, anglebN = 0;
	calvecangle(angleaN, axv, vay);
	calvecangle(anglebN, bxv, vby);
	a.angle = angleaN + alfa, b.angle = anglebN + alfa;
	a.v = sqrt(axv*axv + vay*vay), b.v = sqrt(bxv*bxv + vby*vby);
	skip:
	double morea = 0.0, moreb = 0.0;
	double angle1 = alfa - PI;
	//把重合的长度分成两部分，加到他们的x，y方向累积量上，确保分离
	if (a.v - 0.0 <= Exp&&b.v - 0.0 <= Exp)
	{
		morea = b.mass / part2*more;
		moreb = more - morea;
	}
	else
	{
		double vh = a.v + b.v;//
		morea = a.v / vh*more;
		moreb = more - morea;
	}
	double morexa = morea*cos(angle1);
	double moreya = morea*sin(angle1);
	double morexb = moreb*cos(alfa);
	double moreyb = moreb*sin(alfa);
	//再次确保分离
	if (fabs(morexa - 0.0) >= Exp) morexa += (morexa < 0.0 ? -1.0 : 1.0);
	if (fabs(moreya - 0.0) >= Exp) moreya += (moreya < 0.0 ? -1.0 : 1.0);
	if (fabs(morexb - 0.0) >= Exp) morexb += (morexb < 0.0 ? -1.0 : 1.0);
	if (fabs(moreyb - 0.0) >= Exp) moreyb += (moreyb < 0.0 ? -1.0 : 1.0);
	a.morex = morexa;
	a.morey = moreya;
	b.morex = morexb;
	b.morey = moreyb;
}

void cal(ball *b, int n)//
{
	double dis = 0;
	for (int i = 0; i < n - 1; ++i)
	{
		for (int j = i + 1; j < n; ++j)
		{
			if (b[i].isexist&&b[j].isexist)
			{
				double ugTT = ug*Time*Time;
				double R = b[i].r + b[j].r;
				dis = sqrt((b[i].x - b[j].x)*(b[i].x - b[j].x) + (b[i].y - b[j].y)*(b[i].y - b[j].y)) - R;
				double dis1 = 0.0, dis2 = 0.0;
				if (b[i].move) dis1 = (b[i].v > ug*Time ? b[i].v*Time - 0.5*ugTT : b[i].v*b[i].v / (ug + ug));//
				if (b[j].move) dis2 = (b[j].v > ug*Time ? b[j].v*Time - 0.5*ugTT : b[j].v*b[j].v / (ug + ug));//
				if (dis1 + dis2 > dis)
				{
					double sTime = Time / 10;
					double distance = 0;
					for (int k = 0; k < 10; ++k)
					{
						move(b[i], sTime);
						move(b[j], sTime);
						distance = sqrt((b[i].x - b[j].x)*(b[i].x - b[j].x) + (b[i].y - b[j].y)*(b[i].y - b[j].y));
						if (distance < R)
						{
							//cout << i << " " << j << endl;
							rebound(b[i], b[j], R - distance);
							b[i].move = b[j].move = 1;
						}
					}
					b[i].movech = b[j].movech = 0;
				}
			}
		}
	}
}

hole::hole()
{
	r = 20;
	repaint = 0;
}

bool inhole(ball *b, int ballnum, hole *h, int holenum)
{
	double ugTT = ug*Time*Time;
	for (int i = 0; i < ballnum; ++i)
	{
		if (b[i].isexist&&b[i].r<=20)//hole_r
		{
			double dis = 0.0;
			double dis1 = dis;
			bool isbreak = false;
			for (int j = 0; j < holenum; ++j)
			{
				dis = sqrt((b[i].x - h[j].x)*(b[i].x - h[j].x) + (b[i].y - h[j].y)*(b[i].y - h[j].y));
				dis1 = b[i].v > ug*Time ? b[i].v*Time - 0.5*ugTT : b[i].v*b[i].v / (ug + ug);
				double R = b[i].r + h[j].r;
				dis -= h[j].r;
				if (dis1 > dis)
				{
					double sTime = Time / 10;
					for (int k = 0; k < 10; ++k)
					{
						move(b[i], sTime);
						dis = sqrt((b[i].x - h[j].x)*(b[i].x - h[j].x) + (b[i].y - h[j].y)*(b[i].y - h[j].y));
						if (dis < h[j].r)
						{
							setfillcolor(bkcolor);
							solidcircle(b[i].xl, b[i].yl, b[i].r);
							b[i].isexist = false;
							if (i == mainball) return false;
							++inballnum;
							isbreak = true;
							break;
						}
					}
					b[i].movech = 0;
					if(isbreak) break;
				}
			}
		}
	}
	return true;
}

void repaint(hole *h, int holenum)
{
	setfillcolor(holecolor);
	for (int i = 0; i < holenum; ++i)
	{
		solidcircle(h[i].x, h[i].y, h[i].r);
	}
}

void mouseevent(ball &b)
{
	if (_kbhit())
	{
		char ch = _getch();
		if (ch == ' '&&mlx!=-10)
		{
			setlinecolor(bkcolor);
			line(b.x, b.y, mlx, mly);
			nocancel = false;
		}
	}
	m = GetMouseMsg();
	double vx = m.x - b.x, vy = m.y - b.y;
	double angle = 0.0;
	double dis = sqrt(vx*vx + vy*vy);
	if (m.mkLButton&&dis<=maxlen)
	{
		if (mlx != -10)
		{
			setlinecolor(bkcolor);
			line(b.x, b.y, mlx,mly);
		}
		if (nocancel)
		{
			setlinecolor(linecolor);
			line(b.x, b.y, m.x, m.y);
			mlx = m.x, mly = m.y;
		}
	}
	else if (m.uMsg == WM_LBUTTONUP&&mlx!=-10)
	{
		if (nocancel)
		{
			calvecangle(angle, mlx - b.x, mly - b.y);
			b.angle = angle, b.v = sqrt((mlx - b.x)*(mlx - b.x) + (mly - b.y)*(mly - b.y)) / rate;
			b.move = true;
			setlinecolor(bkcolor);
			line(b.x, b.y, mlx, mly);
			mlx = mly = -10;
		}
		else nocancel = true;
	}
	FlushMouseMsgBuffer();
}

bool judgepair(Pair &p, ball &b)
{
	while (b.angle > PI) b.angle -= (PI + PI);
	while (b.angle <= -PI) b.angle += (PI + PI);
	double angle1, angle2;
	calvecangle(angle1, b.x - p.x[0], b.y - p.y[0]);
	calvecangle(angle2, b.x - p.x[1], b.y - p.y[1]);
	double angle3 = angle2 - angle1;
	while (angle3 < 0) angle3 += (PI + PI);
	while (angle3 >= (PI + PI)) angle3 -= (PI + PI);
	if (angle3 > PI) angle3 = (PI + PI) - angle3;
	double angle4 = angle1 - b.angle, angle5 = angle2 - b.angle;
	while (angle4 < 0) angle4 += (PI + PI);
	while (angle4 >= (PI + PI)) angle4 -= (PI + PI);
	if (angle4 > PI) angle4 = (PI + PI) - angle4;
	while (angle5 < 0) angle5 += (PI + PI);
	while (angle5 >= (PI + PI)) angle5 -= (PI + PI);
	if (angle5 > PI) angle5 = (PI + PI) - angle5;
	if (fabs(angle4 + angle5 - angle3) < Exp*1e3)
	{
		return true;
	}
	return false;
}

bool judgepair1(Pair &p, ball &b)//bug
{
	double angle1, angle2;
	calvecangle(angle1, p.x[0] - b.x, p.y[0] - b.y);
	calvecangle(angle2, p.x[1] - b.x, p.y[1] - b.y);
	if (cos(angle1 - angle2) > 0) return false;
	return true;
}

bool rebound1(cirpoint &s, ball &b)
{
	if (((b.x - s.x)*(b.x - s.x) + (b.y - s.y)*(b.y - s.y)) < b.r*b.r)
	{
		int vecx = b.x - s.x, vecy = b.y - s.y;
		double alfa = 0;
		calvecangle(alfa, vecx, vecy);
		double afterangle = b.angle - alfa;
		if (cos(afterangle) <= 0) b.angle = ((sin(afterangle) > 0) ? PI : -PI) - afterangle + alfa;
		return true;
	}
	return false;
}

void rebound2(int zhangaimark, ball &b)
{
	int k = -zhangaimark;
	int paircnt = 0, pairmark = -1;
	for (int i = 0; i < zhang[k].pairnum; ++i)
	{
		if (zhang[k].p[i].x[0] == zhang[k].p[i].x[1])
		{
			if (judgepair(zhang[k].p[i], b))
			{
				int in = abs(b.x - zhang[k].p[i].x[0]) + b.r;
				if (b.xl <= zhang[k].p[i].x[0]) b.x = b.x - in - in;
				else b.x = b.x + in + in;
				b.angle = PI - b.angle;
				break;
			}
		}
		else
		{
			if (judgepair(zhang[k].p[i], b))
			{
				int in = abs(b.y - zhang[k].p[i].y[0]) + b.r;
				if (b.yl <= zhang[k].p[i].y[0]) b.y = b.y - in - in;
				else b.y = b.y + in + in;
				b.angle =-b.angle;
				break;
			}
		}
	}
}

void iszhangai(zhangai* z, int zhangainum, ball* b, int ballnum)
{
	double sTime = Time / 10;
	for (int i = 0; i < ballnum-solidball; ++i)
	{
		int yes = 1;
		if (b[i].isexist&&b[i].x >= 0 && b[i].y >= 0&&(occupy[b[i].x][b[i].y] == -100|| occupy[b[i].x][b[i].y]>0))//?
		{
			for (int j = 0; j < 10; ++j)
			{
				move(b[i], sTime);
				if (yes&&b[i].x >= 0 && b[i].y >= 0 && b[i].x <= kuan&&b[i].y <= gao&&occupy[b[i].x][b[i].y] > 0)
				{
					if (rebound1(cirp[occupy[b[i].x][b[i].y]], b[i])) yes = 0;
				}
				for (int k = 1; yes&&(k <= zhangainum); ++k)//
				{
					for (int l = 0; l < z[k].pairnum; ++l)
					{
						if (z[k].p[l].x[0] == z[k].p[l].x[1])
						{
							if (abs(b[i].x - z[k].p[l].x[0]) < b[i].r&&b[i].y<=z[k].p[l].y[1]&&b[i].y>=z[k].p[l].y[0])
							{
								int in = b[i].r - abs(b[i].x - z[k].p[l].x[0]);
								if (b[i].x <= z[k].p[l].x[0]) b[i].x = b[i].x - in - in;
								else b[i].x = b[i].x + in + in;
								b[i].angle = PI - b[i].angle;
								yes = 0;
								break;
							}
						}
						else
						{
							if (abs(b[i].y - z[k].p[l].y[0]) < b[i].r&&b[i].x<= z[k].p[l].x[1]&&b[i].x>= z[k].p[l].x[0])
							{
								int in = b[i].r - abs(b[i].y - z[k].p[l].y[0]);
								if (b[i].y <= z[k].p[l].y[0]) b[i].y = b[i].y - in - in;
								else b[i].y = b[i].y + in + in;
								b[i].angle = -b[i].angle;
								yes = 0;
								break;
							}
						}
					}
				}
			}
			b[i].movech = 0;
		}
		else if (b[i].x >= 0 && b[i].y >= 0 && b[i].x <= kuan&&b[i].y <= gao&&b[i].isexist&&occupy[b[i].x][b[i].y] < 0)
		{
			rebound2(occupy[b[i].x][b[i].y], b[i]);
		}
	}
}

void drawer0()
{
	setfillcolor(zhangaicolor);
	solidrectangle(0, 225, 415, 225 + houdu);//左p[0],p[
	solidrectangle(585, 225, kuan, 225 + houdu);//右
	solidrectangle(475, 0, 475 + houdu, 165);//上
	solidrectangle(475, 335, 475 + houdu, gao);//下
}

void inibuilder(int *x, int *y, int zhang_num)
{
	int keycnt = 1;
	memset(occupy, 0, sizeof(occupy));
	setfillcolor(zhangaicolor);
	int firstp = 0;
	int firstpp = 0;
	for (int i = 0; i < zhang_num; ++i)
	{
		solidrectangle(x[firstpp], y[firstpp], x[firstpp + 2], y[firstpp + 2]);
		firstpp += 4;
	}
	for (int i = 1; i <= zhang_num; ++i)
	{
		int jjbegin = (y[firstp] >= 30) ? y[firstp] - 30 : 0;
		int jjend = (y[firstp + 2] <= gao - 30) ? y[firstp + 2] + 30 : gao;
		int kkbegin = (x[firstp] >= 30) ? x[firstp] - 30 : 0;
		int kkend = (x[firstp + 2] <= kuan - 30) ? x[firstp + 2] + 30 : kuan;
		for (int j = jjbegin; j <= jjend; ++j)
		{
			for (int k = kkbegin; k <= kkend; ++k) occupy[k][j] = -100;
		}
		int pair_num = 0;
		if (y[firstp] != 0)
		{
			zhang[i].p[pair_num].x[0] = x[firstp], zhang[i].p[pair_num].x[1] = x[firstp + 1];
			zhang[i].p[pair_num].y[0] = y[firstp], zhang[i].p[pair_num].y[1] = y[firstp + 1];
			pair_num++;
		}
		if (x[firstp + 1] != kuan)
		{
			zhang[i].p[pair_num].x[0] = x[firstp + 1], zhang[i].p[pair_num].x[1] = x[firstp + 2];
			zhang[i].p[pair_num].y[0] = y[firstp + 1], zhang[i].p[pair_num].y[1] = y[firstp + 2];
			pair_num++;
		}
		if (y[firstp + 2] != gao)
		{
			zhang[i].p[pair_num].x[0] = x[firstp + 3], zhang[i].p[pair_num].x[1] = x[firstp + 2];
			zhang[i].p[pair_num].y[0] = y[firstp + 3], zhang[i].p[pair_num].y[1] = y[firstp + 2];
			pair_num++;
		}
		if (x[firstp] != 0)
		{
			zhang[i].p[pair_num].x[0] = x[firstp], zhang[i].p[pair_num].x[1] = x[firstp + 3];
			zhang[i].p[pair_num].y[0] = y[firstp], zhang[i].p[pair_num].y[1] = y[firstp + 3];
			pair_num++;
		}
		zhang[i].pairnum = pair_num;
		firstp += 4;
	}
	firstpp = 0;
	for (int i = 0; i < zhang_num; ++i)
	{
		for (int j = firstpp; j < firstpp + 4; ++j)
		{
			if (x[j] != 0 && x[j] != kuan&&y[j] != 0 && y[j] != gao)
			{
				int tempcnt = 0;
				int coloxl = getpixel(x[j] - 1, y[j]);
				int coloxr = getpixel(x[j] + 1, y[j]);
				int coloys = getpixel(x[j], y[j] - 1);
				int coloyx = getpixel(x[j], y[j] + 1);
				if (coloxl == zhangaicolor) ++tempcnt;
				if (coloxr == zhangaicolor) ++tempcnt;
				if (coloys == zhangaicolor) ++tempcnt;
				if (coloyx == zhangaicolor) ++tempcnt;
				if (tempcnt == 2)
				{
					int kbegin = (coloys == zhangaicolor) ? y[j] : y[j] - 25;//
					int kend = (coloys == zhangaicolor) ? y[j] + 25 : y[j];
					int lbegin = (coloxl == zhangaicolor) ? x[j] : x[j] -25;
					int lend = (coloxl == zhangaicolor) ? x[j] + 25 : x[j];
					cirp[keycnt].x = x[j], cirp[keycnt].y = y[j];
					for (int k = kbegin; k <= kend; ++k)
					{
						for (int l = lbegin; l <= lend; ++l) occupy[l][k] = (((l - x[j])*(l - x[j]) + (k - y[j])*(k - y[j])) <= 625) ? keycnt : -100;
					}
					keycnt++;
				}
			}
		}
		firstpp += 4;
	}
	firstp = 0;
	for (int i = 1; i <= zhang_num; ++i)
	{
		for (int j = y[firstp] + 1; j < y[firstp + 2]; ++j)
			for (int k = x[firstp] + 1; k < x[firstp + 2]; ++k) occupy[k][j] = -i;
		firstp += 4;
	}
}

void builder0()
{
	int x[20] = { 0,415,415,0,475,525,525,475,585,kuan,kuan,585,475,525,525,475 };
	int y[20] = { 225,225,275,275,0,0,165,165,225,225,275,275,335,335,gao,gao };
	zhangnum = 4;
	inibuilder(x, y, 4);
	holenum = 1;
	h[0].x = kuan / 2, h[0].y = gao / 2;
	if (ch != 5)
	{
		ballnum = 14;
		//ini...
		for (int i = 0; i < ballnum; ++i)
		{
			test[i].x = 200, test[i].y = 20, test[i].angle = PI / 3, test[i].isexist = true;
			test[i].morex = test[i].morey = 0;
			test[i].movech = 1, test[i].move = 1, test[i].v = 0;
		}
		bigball = 4;
		solidball = 0;
		for (int i = 10; i < ballnum + bigball; ++i) test[i].mass = 20, test[i].r = 21;
		test[10].x = 225, test[10].y = 112, test[11].x = 725, test[11].y = 112;
		test[12].x = 225, test[12].y = 362, test[13].x = 725, test[13].y = 362;
	}
}

void drawer1()
{
	setfillcolor(zhangaicolor);
	solidrectangle(100, 100, 180, 180);
	solidrectangle(460, 100, 540, 180);
	solidrectangle(820, 100, 900, 180);
	solidrectangle(100, 320, 180, 400);
	solidrectangle(460, 320, 540, 400);
	solidrectangle(820, 320, 900, 400);
}

void builder1()
{
	int x[30] = {100,180,180,100,460,540,540,460,820,900,900,820,100,180,180,100,460,540,540,460,820,900,900,820};
	int y[30] = {100,100,180,180,100,100,180,180,100,100,180,180,320,320,400,400,320,320,400,400,320,320,400,400};
	zhangnum = 6;
	inibuilder(x, y, 6);
	holenum = 2;
	h[0].x = 320, h[0].y = 250, h[1].x = 680, h[1].y = 250;
	if (ch != 5)
	{
		//ini...
		for (int i = 0; i < ballnum; ++i)
		{
			test[i].x = 200, test[i].y = 20, test[i].angle = PI / 3, test[i].isexist = true;
			test[i].morex = test[i].morey = 0;
			test[i].movech = 1, test[i].move = 1, test[i].v = 0;
		}
		bigball = 0;
		solidball = 0;
	}
}

void drawer2()
{
	setfillcolor(zhangaicolor);
	solidrectangle(100, 100, 160, 160);
	solidrectangle(160, 160, 220, 220);
	solidrectangle(220, 220, 280, 280);
	solidrectangle(280, 280, 340, 340);
	solidrectangle(340, 340, 400, 400);
	solidrectangle(600, 100, 660, 160);
	solidrectangle(660, 160, 720, 220);
	solidrectangle(720, 220, 780, 280);
	solidrectangle(780, 280, 840, 340);
	solidrectangle(840, 340, 900, 400);
}

void builder2()
{
	int x[50] = { 100,160,160,100,160,220,220,160,220,280,280,220,280,340,340,280,340,400,400,340,600,660,660,600,660,720,720,660,720,780,780,720,780,840,840,780,840,900,900,840 };
	int y[50] = { 100,100,160,160,160,160,220,220,220,220,280,280,280,280,340,340,340,340,400,400,100,100,160,160,160,160,220,220,220,220,280,280,280,280,340,340,340,340,400,400 };
	zhangnum = 10;
	inibuilder(x, y, zhangnum);
	holenum = 1;
	h[0].x = 500, h[0].y = 250;
	if (ch != 5)
	{
		//ini...
		for (int i = 0; i < ballnum; ++i)
		{
			test[i].x = 200, test[i].y = 20, test[i].angle = PI / 3, test[i].isexist = true;
			test[i].morex = test[i].morey = 0;
			test[i].movech = 1, test[i].move = 1, test[i].v = 0;
		}
		bigball = 0;
		solidball = 0;
	}
}

void drawer3()
{
	setfillcolor(zhangaicolor);
	solidrectangle(175, 120, 825, 220);
	solidrectangle(175, 280, 825, 380);
}

void builder3()
{
	int x[20] = { 175,825,825,175,175,825,825,175 };
	int y[20] = { 120,120,220,220,280,280,380,380 };
	zhangnum = 2;
	inibuilder(x, y, zhangnum);
	holenum = 1;
	h[0].x = 500, h[0].y = 250;
	if (ch != 5)
	{
		ballnum = 12;
		//ini...
		for (int i = 0; i < ballnum; ++i)
		{
			test[i].x = 200, test[i].y = 20, test[i].angle = PI / 3, test[i].isexist = true;
			test[i].morex = test[i].morey = 0;
			test[i].movech = 1, test[i].move = 1, test[i].v = 0;
		}
		bigball = 0;
		solidball = 2;
		test[10].mass = test[11].mass = 1 << 31;
		test[10].x = 70, test[10].y = 250, test[11].x = 930, test[11].y = 250, test[10].r = test[11].r = 30;
	}
}

void buildern()
{
	memset(occupy, 0, sizeof(occupy));
	zhangnum = 0;
	holenum = 0;
	//ini...
	for (int i = 0; i < ballnum; ++i)
	{
		test[i].x = 200, test[i].y = 20, test[i].angle = PI / 3, test[i].isexist = true;
		test[i].morex = test[i].morey = 0;
		test[i].movech = 1, test[i].move = 1, test[i].v = 0;
	}
	bigball = 0;
	solidball = 0;
}

void drawern()
{
	;
}

void choice5()
{
	int ch1;
AGain:
	cout << "选择关卡：1. 2. 3. 4. 5（空）" << endl;
	cin >> ch1;
	if (ch1 < 1 || ch1>5) goto AGain;
	switch (ch1)
	{
	case 1:
	{
		drawer = drawer0;
		builder = builder0;
	}
	break;
	case 2:
	{
		drawer = drawer1;
		builder = builder1;
	}
	break;
	case 3:
	{
		drawer = drawer2;
		builder = builder2;
	}
	break;
	case 4:
	{
		drawer = drawer3;
		builder = builder3;
	}
	break;
	case 5:
	{
		self = 1;
		drawer = drawern;
		builder = buildern;
	}
	break;
	}
Again:
	int ch2;
	cout << "输入球数：";
	cin >> ch2;
	if (ch2 > 499 || ch2 < 1) goto Again;
	ballnum = ch2;
	solidball = bigball = 0;
	//ini
	for (int i = 0; i < ballnum; ++i)
	{
		test[i].mass = 5, test[i].r = 11;
		test[i].x = 200, test[i].y = 20, test[i].angle = PI / 3, test[i].isexist = true;
		test[i].morex = test[i].morey = 0;
		test[i].movech = 1, test[i].move = 1, test[i].v = 0;
	}
	for (int i = 1; i < ballnum; ++i) test[i].x = 500;
	if (ch1 == 1) for (int i = 1; i < ballnum; ++i) test[i].y = 200;
}

void choice()
{
	again:
	cout << "选择关卡：1. 2. 3. 4. 5（bt）" << endl;
	cin >> ch;
	if (ch < 1 || ch>5) goto again;
	ballnum = 10;
	switch (ch)
	{
	case 1:
	{
		drawer = drawer0;
		builder = builder0;
	}
	break;
	case 2:
	{
		drawer = drawer1;
		builder = builder1;
	}
	break;
	case 3:
	{
		drawer = drawer2;
		builder = builder2;
	}
	break;
	case 4:
	{
		drawer = drawer3;
		builder = builder3;
	}
	break;
	case 5:
	{
		choice5();
	}
	break;
	}
}