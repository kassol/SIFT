#include "match.h"
#include <iostream>
using namespace std;


int main(int argc, char* argv[])
{
	match m("D:\\新建文件夹 (2)\\match images\\110270110034.tif", "D:\\新建文件夹 (2)\\match images\\110270111031.tif");
	std::vector<SamePoint> result;
	m.domatch(result);
	system("pause");
	return 0;
}
