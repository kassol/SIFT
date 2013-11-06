#include "match.h"
#include <iostream>
using namespace std;


int main(int argc, char* argv[])
{
	match m("D:\\02-164_50mic.tif", "D:\\02-165_50mic.tif");
	std::vector<SamePoint> result;
	m.domatch(result);
	system("pause");
	return 0;
}
