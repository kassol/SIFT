#ifndef _SIFT_H
#define  _SIFT_H

#import "ImageX.tlb" no_namespace named_guids

#include <list>
#include <vector>

typedef unsigned char uchar;

typedef double pixel_t;



#define PI   3.1415926535897932384626433832795

#define INIT_SIGMA 0.5

#define SIGMA 1.6

#define INTERVALS 3

//比较时的半径，值越小最后取得的特征点越少
#define RATIO 10

#define MAX_INTERPOLATION_STEPS 5

//|D(x^)| < 0.03   0.03极值点太多，增加会使特征点减少
#define DXTHRESHOLD 0.03

#define ORI_HIST_BINS 36

#define ORI_SIGMA_TIMES 1.5

#define ORI_WINDOW_RADIUS 3.0 * ORI_SIGMA_TIMES

#define ORI_SMOOTH_TIMES 2

#define ORI_PEAK_RATIO 0.8

#define FEATURE_ELEMENT_LENGTH 128

#define DESCR_HIST_BINS 8

#define IMG_BORDER 5

#define DESCR_WINDOW_WIDTH 4

#define DESCR_SCALE_ADJUST 3

#define DESCR_MAG_THR 0.2

#define INT_DESCR_FCTR 512.0

#define At(index, x, y) (PyrAt(dog_pyr, (index), (x), (y)))

#define Hat(i, j) (*(H+(i)*3 + (j)))

#define HIat(i, j) (*(H_inve+(i)*3 + (j)))

#define DAt(x, y) (*(data+(y)*step+(x)))

//抛物插值
#define Parabola_Interpolate(l, c, r) (0.5*((l)-(r))/((l)-2.0*(c)+(r)))

#define MINPQ_INIT_NALLOCD 512

#define KDTREE_BBF_MAX_NN_CHKS 200

#define NN_SQ_DIST_RATIO_THR 0.36




//sift算法操作的基本单位-影像块
struct Block{
	pixel_t* pBuffer;
	int nXOrigin;
	int nYOrigin;
	int nXSize;
	int nYSize;
	Block(){pBuffer = NULL; nXOrigin = 0; nYOrigin = 0; nXSize = 0; nYSize = 0;}
	Block(pixel_t* buffer, int ox, int oy, int nx, int ny){pBuffer = buffer; nXOrigin = ox; nYOrigin = oy; nXSize = nx; nYSize = ny;}
};



struct Keypoint
{
	int octave; //关键点所在组
	int interval;// 关键点所在层

	double offset_interval;//调整后的层的增量

	int x; //x,y坐标,根据octave和interval可取的层内图像
	int y;

	//scale = sigma0*pow(2.0, o+s/S)
	double scale; //空间尺度坐标

	double dx; //特征点坐标，该坐标被缩放成原图像大小 
	double dy;

	double offset_x;
	double offset_y;

	//高斯金字塔组内各层尺度坐标，不同组的相同层的sigma值相同
	//关键点所在组的组内尺度
	double octave_scale; //offset_i;

	double ori;//方向

	int descr_length;
	double descriptor[FEATURE_ELEMENT_LENGTH]; //描述子

	double val;//极值
	void* feature_data;
};

//提取特征点
void sift(pixel_t* init_pbuf, std::vector<Keypoint>&features,int nxsize, int nysize);


//k-d树的节点
struct kd_node
{
	int ki;										//描述子中关键位的索引
	double kv;								// 描述子中关键位的值
	int leaf;									//如果是叶节点为1，否则为0
	Keypoint* features;				//该节点的特征点
	int n;										//该节点特征点的个数
	kd_node* kd_left;					//左孩子  
	kd_node* kd_right;				//右孩子
};

struct bbf_data
{
	double d;
	void* old_data;
};

struct pq_node
{
	void* data;
	int key;
};


struct min_pq
{
	pq_node* pq_array;		//优先队列
	int nallocd;					//分配元素的个数
	int n;								//队列中元素的个数
};




//k-d树的基本操作
kd_node* kdtree_build(Keypoint* features, int n);

void kdtree_release(kd_node* kd_root);

int kdtree_bbf_knn(kd_node* kd_root, Keypoint* feat, int k, Keypoint*** nbrs, int max_nn_chks);

double descr_dist_sq(Keypoint* f1, Keypoint* f2);



//影像中的同名点
struct SamePoint
{
	double lx;
	double ly;
	double rx;
	double ry;

	SamePoint(){lx = 0; ly = 0; rx = 0; ry = 0;}
	SamePoint(double mlx, double mly, double mrx, double mry)
	{lx = mlx; ly = mly; rx = mrx; ry = mry;}
};


class Sift{
public:
	Sift()
	{
		CoInitialize(NULL);
		hthread = NULL;
		m_pImage = NULL;
		m_bFinish = false;
		::CoCreateInstance(CLSID_ImageDriver, NULL, CLSCTX_ALL, IID_IImage, (void**)&m_pImage);
	}
	~Sift()
	{
		m_pImage->Release();
		if (hthread != NULL)
		{
			delete []hthread;
			hthread = NULL;
		}
		::CoUninitialize();
	}
	void fetchFeatures(const char* szPathName);
	Keypoint* getFeatures();

public:
	std::list<Block> m_listBlock;
	std::list<Keypoint> m_listKeyPoint;
	bool m_bFinish;
	HANDLE hmutex;
	HANDLE hmutex2;

private:
	IImage* m_pImage;
	HANDLE* hthread;
};


#endif
