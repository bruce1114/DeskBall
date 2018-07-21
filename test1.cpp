#define _CRT_SECURE_NO_WARNINGS
#include<iostream>
#include<cmath>
#include<conio.h>
#include"ball.h"
#include<ctime>
using namespace std;

int main()
{
	cout << "游戏结束时，输入y继续，否则退出" << endl;
	again:
	choice();

	initgraph(kuan, gao+extra);
	setfillcolor(bkcolor);
	solidrectangle(-1, -1, kuan, gao);

	builder();
	drawer();
	test[0].y = 20;
	test[0].x = 20;

	setfillcolor(mainbcolor);
	solidcircle(test[0].x, test[0].y, test[0].r), test[0].xl = test[0].x, test[0].yl = test[0].y;

	setfillcolor(otherbcolor);
	for (int i = mainball + 1; i < ballnum; ++i) solidcircle(test[i].x, test[i].y, test[i].r), test[i].xl = test[i].x, test[i].yl = test[i].y;
	setfillcolor(zhangaicolor);
	for (int i = ballnum - solidball; i < ballnum; ++i) solidcircle(test[i].x, test[i].y, test[i].r);

	setfillcolor(holecolor);
	for (int i = 0; i < holenum; ++i) solidcircle(h[i].x, h[i].y, h[i].r);

	BeginBatchDraw();
	while (true)
	{
		Sleep(sleeptime);
		setfillcolor(bkcolor);
		for(int i=0;i<ballnum;++i) if(test[i].isexist) solidcircle(test[i].xl, test[i].yl, test[i].r);

		drawer();
		//repaint
		setfillcolor(holecolor);
		for (int i = 0; i < holenum; ++i) solidcircle(h[i].x, h[i].y, h[i].r);

		cal(test, ballnum);

		iszhangai(zhang, zhangnum, test, ballnum);

		for(int i=0;i<ballnum-solidball;++i) move(test[i],Time);
		for (int i = ballnum - solidball; i < ballnum; ++i) test[i].x = test[i].xl, test[i].y = test[i].yl;

		setfillcolor(mainbcolor);
		if(test[mainball].isexist) solidcircle(test[0].x, test[0].y, test[0].r), test[0].xl = test[0].x, test[0].yl = test[0].y;

		setfillcolor(otherbcolor);
		for (int i = mainball + 1; i < ballnum-solidball; ++i) if(test[i].isexist) solidcircle(test[i].x, test[i].y, test[i].r), test[i].xl = test[i].x, test[i].yl = test[i].y;
		setfillcolor(zhangaicolor);
		for (int i = ballnum - solidball; i < ballnum; ++i) solidcircle(test[i].x, test[i].y, test[i].r);

		if (!test[0].move&&test[0].isexist&&MouseHit()) mouseevent(test[0]);
		FlushBatchDraw();
		if (!inhole(test, ballnum, h, holenum))
		{
			//repaint
			setfillcolor(holecolor);
			for (int i = 0; i < holenum; ++i) solidcircle(h[i].x, h[i].y, h[i].r);
			FlushBatchDraw();
			break;
		}
		if (inballnum == ballnum - 1-bigball-solidball)
		{
			//repaint
			setfillcolor(holecolor);
			for (int i = 0; i < holenum; ++i) solidcircle(h[i].x, h[i].y, h[i].r);
			FlushBatchDraw();
			break;
		}
		if (self) test[0].v = 5, test[0].move = 1;
		if (_kbhit())
		{
			_getch();

		}
	}
	inballnum = 0;
	char ch;
	ch = _getch();
	closegraph();
	self = 0;
	if (ch == 'y') goto again;
	return 0;
}