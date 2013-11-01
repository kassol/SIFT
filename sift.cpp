#include "sift.h"
#include <Windows.h>
#include <algorithm>
#include <time.h>

const int nBlockSize = 512;


//向上采样
void UpSample(const pixel_t* src, pixel_t** dst, int &nxsize, int &nysize);

//向下采样
void DownSample(const pixel_t* src, pixel_t**dst, int &nxsize, int &nysize);

//高斯模糊
void GaussianSmooth(const pixel_t* src, pixel_t** dst, int nxsize, int nysize, double sigma);

//建立高斯金字塔
void GaussianPyramid(pixel_t* src, std::vector<Block>&gauss_pyr, int nxsize, int nysize, int octaves, 
         int intervals = INTERVALS, double sigma = SIGMA);

//建立差分金字塔
void DogPyramid(const std::vector<Block>&gauss_pyr, std::vector<Block>&dog_pyr, 
         int octaves, int intervals = INTERVALS);

void Sub(const pixel_t* a, const pixel_t* b, pixel_t** c, int nxsize, int nysize);

//提取极值点
void FetchExtrema(const std::vector<Block>&dog_pyr, std::vector<Keypoint>&extrema,
	     int octaves, int intervals = INTERVALS);

bool IsExtremum(int x, int y, const std::vector<Block>& dog_pyr, int index);

//去除边缘点和不稳定点
Keypoint* InterploationExtremum(int x, int y, const std::vector<Block>& dog_pyr, int index, int octave, 
                   int interval, double dxthreshold = DXTHRESHOLD);
double GetFabsDx(int x, int y, const std::vector<Block>& dog_pyr, int index, const double* offset_x);
void GetOffsetX(int x, int y, const std::vector<Block>& dog_pyr, int index, double *offset_x);
void Hessian3D(int x, int y, const std::vector<Block>& dog_pyr, int index, double *H);
double PyrAt(const std::vector<Block>& pyr, int index, int x, int y);
bool Inverse3D(const double *H, double *H_inve);
void DerivativeOf3D(int x, int y, const std::vector<Block>& dog_pyr, int index, double *dx);
bool passEdgeResponse(int x, int y, const std::vector<Block>& dog_pyr, int index, double r = RATIO);

//计算尺度
void CalculateScale(std::vector<Keypoint>& features, double sigma = SIGMA, 
	     int intervals = INTERVALS);
void HalfFeatures(std::vector<Keypoint>& features);

//分配方向
void OrientationAssignment(std::vector<Keypoint>& extrema, std::vector<Keypoint>& features, 
	     const std::vector<Block>& gauss_pyr);

//计算方向直方图
double* CalculateOrientationHistogram(const Block& gauss, int x, int y, int bins, int radius, double sigma);
bool CalcGradMagOri(const Block& gauss, int x, int y, double& mag, double& ori);
void GaussSmoothOriHist(double *hist, int n);
double DominantDirection(double *hist, int n);
void CalcOriFeatures(const Keypoint& keypoint, std::vector<Keypoint>& features, 
	     const double *hist, int n, double mag_thr);
void CopyKeypoint(const Keypoint& src, Keypoint& dst);
void DescriptorRepresentation(std::vector<Keypoint>& features, 
	     const std::vector<Block>& gauss_pyr, int bins, int width);
double*** CalculateDescrHist(const Block& gauss, int x, int y, 
	                double octave_scale, double ori, int bins, int width);
void InterpHistEntry(double ***hist, double xbin, double ybin, double obin, double mag, int bins, int d);
void HistToDescriptor(double ***hist, int width, int bins, Keypoint& feature);
void NormalizeDescr(Keypoint& feat);
bool FeatureCmp(Keypoint& f1, Keypoint& f2);

//k-d tree、BBF搜索
kd_node* kd_node_init(Keypoint* features, int n);

void expand_kd_node_subtree(kd_node* node);
void assign_part_key(kd_node* kd_node);
double median_select(double* array, int n);
double rank_select(double* array, int n, int r);
void partition_features(kd_node* node);
void insertion_sort(double* array, int n);
int partition_array(double* array, int n, double pivot);
min_pq* minpq_init();
int minpq_insert(min_pq* pq, void* data, int key);
void* minpq_extract_min(min_pq* pq);
void minpq_release(min_pq** pq);
kd_node* explore_to_leaf(kd_node* node, Keypoint* feat, min_pq* min_pq);
int insert_into_nbr_array(Keypoint* feat, Keypoint** nbrs, int n, int k);
int array_double(void** array, int n, int size);
void decrease_pq_node_key(pq_node* pq_array, int i, int key);
void restore_minpq_order(pq_node* pq_array, int i, int n);




static inline int parent(int i)
{
	return (i - 1) / 2;
}


static inline int right(int i)
{
	return 2 * i + 2;
}

static inline int left(int i)
{
	return 2 * i + 1;
}

void sift(pixel_t* init_pbuf, std::vector<Keypoint>&features,int nxsize, int nysize)
{
	pixel_t* init_grey = NULL;
	UpSample(init_pbuf, &init_grey, nxsize, nysize);
	delete []init_pbuf;
	init_pbuf = NULL;
	pixel_t* temp;
	double  sigma_init = sqrt(SIGMA * SIGMA - (INIT_SIGMA * 2) * (INIT_SIGMA * 2));
	GaussianSmooth(init_grey, &temp, nxsize, nysize, sigma_init);
	delete []init_grey;
	init_grey = temp;
	temp = NULL;
	int octaves = int(log((double)min(nxsize, nysize))/log(2.0) - 2);
	
	std::vector<Block> gauss_pyr;
	GaussianPyramid(init_grey, gauss_pyr, nxsize, nysize, octaves, INTERVALS, SIGMA);

	std::vector<Block> dog_pyr;
	DogPyramid(gauss_pyr, dog_pyr, octaves);

	std::vector<Keypoint> extrema;
	FetchExtrema(dog_pyr, extrema, octaves);

	CalculateScale(extrema);
	HalfFeatures(extrema);

	OrientationAssignment(extrema, features, gauss_pyr);

	extrema.swap(std::vector<Keypoint>());

	DescriptorRepresentation(features, gauss_pyr, DESCR_HIST_BINS, DESCR_WINDOW_WIDTH);
	
	sort(features.begin(), features.end(), FeatureCmp);
	
	std::vector<Block>::iterator temIte = gauss_pyr.begin();
	while(temIte != gauss_pyr.end())
	{
		delete []temIte->pBuffer;
		++temIte;
	}
	gauss_pyr.clear();
	temIte = dog_pyr.begin();
	while(temIte != dog_pyr.end())
	{
		delete []temIte->pBuffer;
		++temIte;
	}
	dog_pyr.clear();
}


void threadFunc(LPVOID param)
{
	Sift* si = (Sift*)param;
	while(!si->m_listBlock.empty() || si->m_bFinish == false)
	{
		if (si->m_listBlock.empty())
		{
			Sleep(10);
			continue;
		}
		WaitForSingleObject(si->hmutex, INFINITE);
		if (si->m_listBlock.empty())
		{
			ReleaseMutex(si->hmutex);
			Sleep(100);
			continue;
		}
		Block temblock = si->m_listBlock.front();
		si->m_listBlock.pop_front();
		ReleaseMutex(si->hmutex);

		//DOG,feature.
		std::vector<Keypoint> vecKeyPoint;
		sift(temblock.pBuffer, vecKeyPoint, temblock.nXSize, temblock.nYSize);
		while(!vecKeyPoint.empty())
		{
			Keypoint keypoint = vecKeyPoint.back();
			vecKeyPoint.pop_back();
			keypoint.dx += temblock.nXOrigin;
			keypoint.dy += temblock.nYOrigin;
			WaitForSingleObject(si->hmutex2, INFINITE);
			si->m_listKeyPoint.push_back(keypoint);
			ReleaseMutex(si->hmutex2);
		}
		vecKeyPoint.swap(std::vector<Keypoint>());
	}
}

void Sift::fetchFeatures(const char* szPathName)
{
	SYSTEM_INFO sys;
	GetSystemInfo(&sys);
	int n = (int)sys.dwNumberOfProcessors;
	hthread = new HANDLE[n];
	hmutex = CreateMutex(NULL, FALSE, NULL);
	hmutex2 = CreateMutex(NULL, FALSE, NULL);
	for (int i = 0; i < n; ++i)
	{
		hthread[i] = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)threadFunc, LPVOID(this),0, NULL);
	}
	m_pImage->Open(szPathName, modeRead);
	int nXSize, nYSize, nBands;
	m_pImage->GetCols(&nXSize);
	m_pImage->GetRows(&nYSize);
	m_pImage->GetBandNum(&nBands);
	const int nStepSize = 512;
	int nx, ny;
	if (nBands == 1)
	{
		for (int j = 0; j <nYSize; j += nStepSize)
		{
			for (int i = 0; i < nXSize; i += nStepSize)
			{
				nx = nBlockSize;
				ny = nBlockSize;
				if (i+nStepSize >= nXSize)
				{
					nx = nXSize-i;
				}
				if (j+nStepSize >= nYSize)
				{
					ny = nYSize-j;
				}
				uchar* pBuf = new uchar[nx*ny];
				m_pImage->ReadImg(i, j, i+nx, j+ny, pBuf, nx, ny, nBands, 0, 0, nx, ny, -1, 0);
				pixel_t* pBlock = new pixel_t[nx*ny];
				for (int n = 0; n < nx*ny; ++n)
				{
					pBlock[n] = pBuf[n]/255.0;
				}
				delete []pBuf;
				pBuf = NULL;
				while(m_listBlock.size()>10)
				{
					Sleep(100);
				}
				m_listBlock.push_back(Block(pBlock, i, j, nx, ny));
			}
		}
	}
	else if (nBands >= 3)
	{
		for (int j = 0; j <nYSize; j += nStepSize)
		{
			for (int i = 0; i < nXSize; i += nStepSize)
			{
				nx = nBlockSize;
				ny = nBlockSize;
				if (i+nStepSize >= nXSize)
				{
					nx = nXSize-i;
				}
				if (j+nStepSize >= nYSize)
				{
					ny = nYSize-j;
				}
				uchar* pBuf = new uchar[nx*ny*3];
				m_pImage->ReadImg(i, j, i+nx, j+ny, pBuf, nx, ny, 3, 0, 0, nx, ny, 0, 0);
				m_pImage->ReadImg(i, j, i+nx, j+ny, pBuf, nx, ny, 3, 0, 0, nx, ny, 1, 1);
				m_pImage->ReadImg(i, j, i+nx, j+ny, pBuf, nx, ny, 3, 0, 0, nx, ny, 2, 2);
				pixel_t* pBlock = new pixel_t[nx*ny];
				for (int n = 0; n < ny; ++n)
				{
					for (int m = 0; m < nx; ++m)
					{
						pBlock[n*nx+m] =(pBuf[n*nx*3+m*3]+pBuf[n*nx*3+m*3+1]+pBuf[n*nx*3+m*3+2])/3.0/255.0; 
					}
				}
				delete []pBuf;
				pBuf = NULL;
				while(m_listBlock.size()>10)
				{
					Sleep(100);
				}
				m_listBlock.push_back(Block(pBlock, i, j, nx, ny));
			}
		}
	}
	else
	{
		return;
	}
	m_bFinish = true;
	WaitForMultipleObjects(n, hthread, TRUE, INFINITE);
}



void UpSample(const pixel_t* src, pixel_t** dst, int &nxsize, int &nysize)
{
	int nWidth = nxsize*2;
	int nHeight = nysize*2;
	*dst = new pixel_t[nWidth*nHeight];
	
	for (int j = 0; j < nysize-1; ++j)
	{
		for (int i =0; i < nxsize-1; ++i)
		{
			(*dst)[j*2*nWidth+i*2] = src[j*nxsize+i];
			(*dst)[j*2*nWidth+i*2+1] = (src[j*nxsize+i]+src[j*nxsize+i+1])/2.0;
			(*dst)[(j*2+1)*nWidth+i*2] = (src[j*nxsize+i]+src[(j+1)*nxsize+i])/2.0;
			(*dst)[(j*2+1)*nWidth+i*2+1] = (src[j*nxsize+i]+src[j*nxsize+i+1]  \
				+src[(j+1)*nxsize+i]+src[(j+1)*nxsize+i+1])/4.0;
		}
	}
	//最后两行两列
	for (int i = 0; i < nWidth; ++i)
	{
		(*dst)[(nHeight-2)*nWidth+i] = src[(nysize-1)*nxsize+i/2];
		(*dst)[(nHeight-1)*nWidth+i] = src[(nysize-1)*nxsize+i/2];
	}
	for (int i = 0; i < nHeight; ++i)
	{
		(*dst)[i*nWidth+nWidth-2] = src[i/2*nxsize+nxsize-1];
		(*dst)[i*nWidth+nWidth-1] = src[i/2*nxsize+nxsize-1];
	}
	nxsize = nWidth;
	nysize = nHeight;
}

void DownSample(const pixel_t* src, pixel_t**dst, int &nxsize, int &nysize)
{
	int nWidth = nxsize/2;
	int nHeight = nysize/2;
	*dst = new pixel_t[nWidth*nHeight];
	int m = 0, n = 0;

	for (int j = 0; j < nysize; j += 2, ++n)
	{
		m = 0;
		for (int i = 0; i < nxsize; i += 2, ++m)
		{
			if (m < nWidth && n < nHeight)
			{
				(*dst)[n*nWidth+m] = src[j*nxsize+i];
			}
		}
	}

	nxsize = nWidth;
	nysize = nHeight;
}

void GaussianSmooth(const pixel_t* src, pixel_t** dst, int nxsize, int nysize, double sigma)
{
	*dst = new pixel_t[nxsize*nysize];
	pixel_t* temData = new pixel_t[nxsize*nysize];
	sigma = sigma > 0 ? sigma : -sigma;
	int ksize = int(sigma*3+0.5)*2+1;
	double* kernel = new double[ksize];

	double scale = -0.5/(sigma*sigma);
	double cons = 1/sqrt(-scale/PI);

	double sum = 0;
	int kcenter = ksize/2;

	for (int i = 0; i < ksize; ++i)
	{
		int x = i -kcenter;
		kernel[i] = cons*exp(x*x*scale);
		sum += kernel[i];
	}

	for (int i = 0; i < ksize; ++i)
	{
		kernel[i] /= sum;
	}

	for (int y = 0; y < nysize; ++y)
	{
		for (int x = 0; x < nxsize; ++x)
		{
			double mul = 0;
			sum = 0;
			for (int i = -kcenter; i <= kcenter; ++i)
			{
				if (x+i >= 0 && x+i < nxsize)
				{
					mul += src[y*nxsize+x+i]*kernel[kcenter+i];
					sum += kernel[kcenter+i];
				}
			}
			temData[y*nxsize+x] = mul/sum;
		}
	}

	for (int x = 0; x < nxsize; ++x)
	{
		for (int y = 0; y < nysize; ++y)
		{
			double mul = 0;
			sum = 0;
			for(int i = -kcenter; i <= kcenter; ++i)
			{
				if (y+i >= 0 && y +i < nysize)
				{
					mul += temData[(y+i)*nxsize+x]*kernel[kcenter+i];
					sum += kernel[kcenter+i];
				}
			}
			(*dst)[y*nxsize+x] = mul/sum;
		}
	}

	delete []temData;
	temData = NULL;
	delete []kernel;
	kernel = NULL;
}

void GaussianPyramid(pixel_t* src, std::vector<Block>&gauss_pyr, int nxsize, int nysize, int octaves, 
                                         int intervals /* = INTERVALS */, double sigma /* = SIGMA */)
{
	double* sigmas = new double[intervals+3];
	double k = pow(2.0, 1.0/intervals);

	sigmas[0] = sigma;
	double sig_prev, sig_total;

	for(int i = 1; i < intervals + 3; ++i)
	{
		sig_prev = pow(k, i - 1)*sigma;
		sig_total = sig_prev*k;
		sigmas[i] = sqrt(sig_total*sig_total-sig_prev*sig_prev);
	}

	for (int o = 0; o < octaves; ++o)
	{
		for (int i = 0; i < intervals+3; ++i)
		{
			Block block;
			if (o == 0 && i == 0)
			{
				block.nXSize = nxsize;
				block.nYSize = nysize;
				block.pBuffer = src;
			}
			else if (i == 0)
			{
				pixel_t* temp = gauss_pyr[(o-1)*(intervals+3)+intervals].pBuffer;
				int nx = gauss_pyr[(o-1)*(intervals+3)+intervals].nXSize;
				int ny = gauss_pyr[(o-1)*(intervals+3)+intervals].nYSize;
				DownSample(temp, &block.pBuffer, nx, ny);
				block.nXSize = nx;
				block.nYSize = ny;
			}
			else
			{
				pixel_t* temp = gauss_pyr[o*(intervals+3)+i-1].pBuffer;
				int nx = gauss_pyr[o*(intervals+3)+i-1].nXSize;
				int ny = gauss_pyr[o*(intervals+3)+i-1].nYSize;
				GaussianSmooth(temp, &block.pBuffer, nx, ny, sigmas[i]);
				block.nXSize = nx;
				block.nYSize = ny;
			}
			gauss_pyr.push_back(block);
		}
	}

	delete []sigmas;
	sigmas = NULL;
}

void DogPyramid(const std::vector<Block>&gauss_pyr, std::vector<Block>&dog_pyr, 
                                 int octaves, int intervals /* = INTERVALS */)
{
	for (int o = 0; o < octaves; ++o)
	{
		for (int i = 1; i < intervals+3; ++i)
		{
			Block block;
			int nx = gauss_pyr[o*(intervals+3)+i].nXSize;
			int ny = gauss_pyr[o*(intervals+3)+i].nYSize;
			Sub(gauss_pyr[o*(intervals+3)+i].pBuffer, gauss_pyr[o*(intervals+3)+i-1].pBuffer, 
                   &block.pBuffer, nx, ny);
			block.nXSize = nx;
			block.nYSize = ny;
			dog_pyr.push_back(block);
		}
	}
}

void Sub(const pixel_t* a, const pixel_t* b, pixel_t** c, int nxsize, int nysize)
{
	*c = new pixel_t[nxsize*nysize];

	for (int i = 0; i < nxsize*nysize; ++i)
	{
		(*c)[i] = b[i]-a[i];
	}
}

void FetchExtrema(const std::vector<Block>&dog_pyr, std::vector<Keypoint>&extrema, 
                                   int octaves, int intervals /* = INTERVALS */)
{
	long int dd = 0, cc1 = 0, cc2 = 0, cc3 = 0, cc0 = 0, cc00=0;
	double thresh = 0.5*DXTHRESHOLD/intervals;

	for (int o = 0; o < octaves; ++o)
	{
		for (int i = 1; i < intervals+1; ++i)
		{
			int index = o*(intervals+2)+i;
			pixel_t* data = dog_pyr[index].pBuffer;
			int step = dog_pyr[index].nXSize;
			for (int y = IMG_BORDER; y < dog_pyr[index].nYSize-IMG_BORDER; ++y)
			{
				for (int x = IMG_BORDER; x < dog_pyr[index].nXSize-IMG_BORDER; ++x)
				{
					++cc00;
					pixel_t val = data[y*step+x];
					if (fabs(val) > thresh)
					{
						++cc0;
						if (IsExtremum(x, y, dog_pyr, index))
						{
							++cc1;
							Keypoint *extrmum = InterploationExtremum(x, y, dog_pyr, index, o, i);

							if (extrmum)
							{
								++cc2;
								if (passEdgeResponse(extrmum->x, extrmum->y, dog_pyr, index))
								{
									extrmum->val = data[extrmum->y*step+extrmum->x];
									++cc3;
									extrema.push_back(*extrmum);
								}
								delete extrmum;
							}
						}
					}
				}
			}
		}
	}
}

bool IsExtremum(int x, int y, const std::vector<Block>& dog_pyr, int index)
{
	pixel_t* data = dog_pyr[index].pBuffer;
	int step = dog_pyr[index].nXSize;
	pixel_t val = data[y*step+x];

	if (val>0)
	{
		for (int n = -1; n <= 1; ++n)
		{
			for (int j = -1; j <= 1; ++j)
			{
				for (int i = -1; i <= 1; ++i)
				{
					if (val < dog_pyr[index+n].pBuffer[(y+j)*step+x+i])
					{
						return false;
					}
				}
			}
		}
	}
	else
	{
		for (int n = -1; n <= 1; ++n)
		{
			for (int j = -1; j <= 1; ++j)
			{
				for (int i = -1; i <= 1; ++i)
				{
					if (val > dog_pyr[index+n].pBuffer[(y+j)*step+x+i])
					{
						return false;
					}
				}
			}
		}
	}

	return true;
}

Keypoint* InterploationExtremum(int x, int y, const std::vector<Block>& dog_pyr, int index, int octave, 
                                                               int interval, double dxthreshold /* = DXTHRESHOLD */)
{
	double offset_x[3] = {0};

	const Block& block = dog_pyr[index];
	int idx = index;
	int intvl = interval;
	int i = 0;

	while(i < MAX_INTERPOLATION_STEPS)
	{
		GetOffsetX(x, y, dog_pyr, idx, offset_x);
		if(fabs(offset_x[0]) < 0.5 && fabs(offset_x[1]) < 0.5 && fabs(offset_x[2]) < 0.5)
		{
			break;
		}

		x += int(offset_x[0]+0.5);
		y += int(offset_x[1]+0.5);
		interval += int(offset_x[2]+0.5);

		idx = index-intvl+interval;

		if(interval < 1 || interval > INTERVALS ||
			x >= block.nXSize-1 || x < 2 ||
			y >= block.nYSize-1 || y < 2) 
		{
			return NULL;
		}
		++i;
	}

	if(i >= MAX_INTERPOLATION_STEPS)
	{
		return NULL;
	}

	if(GetFabsDx(x, y, dog_pyr, idx, offset_x) < dxthreshold/INTERVALS)
	{
		return NULL;
	}

	Keypoint *keypoint = new Keypoint;

	keypoint->x = x;
	keypoint->y = y;

	keypoint->offset_x = offset_x[0];
	keypoint->offset_y = offset_x[1];

	keypoint->interval = interval;
	keypoint->offset_interval = offset_x[2];

	keypoint->octave = octave;

	keypoint->dx = (x + offset_x[0])*pow(2.0, octave);
	keypoint->dy = (y + offset_x[1])*pow(2.0, octave);

	return keypoint;
}

double GetFabsDx(int x, int y, const std::vector<Block>& dog_pyr, int index, const double* offset_x)
{
	double dx[3];
	DerivativeOf3D(x, y, dog_pyr, index, dx);

	double term = 0.0;

	for (int i = 0; i < 3; ++i)
	{
		term += dx[i]*offset_x[i];
	}

	pixel_t* data = dog_pyr[index].pBuffer;
	int step = dog_pyr[index].nXSize;
	pixel_t val = data[y*step+x];

	return fabs(val+0.5*term);
}

void GetOffsetX(int x, int y, const std::vector<Block>& dog_pyr, int index, double *offset_x)
{
	double H[9], H_inve[9]= {0};
	Hessian3D(x, y, dog_pyr, index, H);
	Inverse3D(H, H_inve);

	double dx[3];
	DerivativeOf3D(x, y, dog_pyr, index, dx);

	for(int i = 0; i < 3; i ++)
	{
		offset_x[i] = 0.0;
		for(int j = 0; j < 3; j++)
		{
			offset_x[i] += H_inve[i*3 + j] * dx[j];
		}
		offset_x[i] = -offset_x[i];
	}
}

void Hessian3D(int x, int y, const std::vector<Block>& dog_pyr, int index, double *H)
{
	double val, Dxx, Dyy, Dss, Dxy, Dxs, Dys;

	val = At(index, x, y);

	Dxx = At(index, x+1, y) + At(index, x-1, y) - 2*val;
	Dyy = At(index, x, y+1) + At(index, x, y-1) - 2*val;
	Dss = At(index+1, x, y) + At(index-1, x, y) - 2*val;

	Dxy = (At(index, x+1, y+1) + At(index, x-1, y-1)
		- At(index, x+1, y-1) - At(index, x-1, y+1))/4.0;

	Dxs = (At(index+1, x+1, y) + At(index-1, x-1, y)
		- At(index-1, x+1, y) - At(index+1, x-1, y))/4.0;

	Dys = (At(index+1, x, y+1) + At(index-1, x, y-1)
		- At(index+1, x, y-1) - At(index-1, x, y+1))/4.0;

	Hat(0, 0) = Dxx;
	Hat(1, 1) = Dyy;
	Hat(2, 2) = Dss;

	Hat(1, 0) = Hat(0, 1) = Dxy;
	Hat(2, 0) = Hat(0 ,2) = Dxs;
	Hat(2, 1) = Hat(1, 2) = Dys;
}

double PyrAt(const std::vector<Block>& pyr, int index, int x, int y)
{
	pixel_t* data = pyr[index].pBuffer;
	int step = pyr[index].nXSize;
	return data[y*step+x];
}

bool Inverse3D(const double *H, double *H_inve)
{
	double A = Hat(0, 0)*Hat(1, 1)*Hat(2, 2) 
		+ Hat(0, 1)*Hat(1, 2)*Hat(2, 0)
		+ Hat(0, 2)*Hat(1, 0)*Hat(2, 1)
		- Hat(0, 0)*Hat(1, 2)*Hat(2, 1)
		- Hat(0, 1)*Hat(1, 0)*Hat(2, 2)
		- Hat(0, 2)*Hat(1, 1)*Hat(2, 0);

	if(fabs(A) < 1e-10)
	{
		return false;
	}

	HIat(0, 0) = Hat(1, 1) * Hat(2, 2) - Hat(2, 1)*Hat(1, 2);
	HIat(0, 1) = -(Hat(0, 1) * Hat(2, 2) - Hat(2, 1) * Hat(0, 2));
	HIat(0, 2) = Hat(0, 1) * Hat(1, 2) - Hat(0, 2)*Hat(1, 1);

	HIat(1, 0) = Hat(1, 2) * Hat(2, 0) - Hat(2, 2)*Hat(1, 0);
	HIat(1, 1) = -(Hat(0, 2) * Hat(2, 0) - Hat(0, 0) * Hat(2, 2));
	HIat(1, 2) = Hat(0, 2) * Hat(1, 0) - Hat(0, 0)*Hat(1, 2);

	HIat(2, 0) = Hat(1, 0) * Hat(2, 1) - Hat(1, 1)*Hat(2, 0);
	HIat(2, 1) = -(Hat(0, 0) * Hat(2, 1) - Hat(0, 1) * Hat(2, 0));
	HIat(2, 2) = Hat(0, 0) * Hat(1, 1) - Hat(0, 1)*Hat(1, 0);

	for(int i = 0; i < 9; i++)
	{
		H_inve[i] /= A;
	}

	return true;
}

void DerivativeOf3D(int x, int y, const std::vector<Block>& dog_pyr, int index, double *dx)
{
	dx[0] = (At(index, x+1, y)-At(index, x-1, y))/2.0;
	dx[1] = (At(index, x, y+1)-At(index, x, y-1))/2.0;
	dx[2] = (At(index+1, x, y)-At(index-1, x, y))/2.0;
}

bool passEdgeResponse(int x, int y, const std::vector<Block>& dog_pyr, int index, 
	                                          double r /* = RATIO */)
{
	pixel_t* data = dog_pyr[index].pBuffer;
	int step = dog_pyr[index].nXSize;
	pixel_t val = data[y*step+x];

	double Dxx, Dyy, Dxy;
	double Tr_h, Det_h;

	Dxx = DAt(x+1, y) + DAt(x-1, y) - 2*val;
	Dyy = DAt(x, y+1) + DAt(x, y-1) - 2*val;
	Dxy = (DAt(x+1, y+1) + DAt(x-1,y-1) - DAt(x-1, y+1) - DAt(x+1, y-1))/4.0;

	Tr_h = Dxx+Dyy;
	Det_h = Dxx*Dyy - Dxy*Dxy;

	if(Det_h <=0)
	{
		return false;
	}

	if(Tr_h*Tr_h/Det_h < (r+1)*(r+1)/r)
	{
		return true;
	}

	return false;
}

void CalculateScale(std::vector<Keypoint>& features, double sigma /* = SIGMA */, 
	                                 int intervals /* = INTERVALS */)
{
	double intvl = 0;

	for (size_t i = 0; i < features.size(); ++i)
	{
		intvl = features[i].interval+features[i].offset_interval;
		features[i].scale = sigma*pow(2.0, features[i].octave+intvl/intervals);
		features[i].octave_scale = sigma*pow(2.0, intvl/intervals);
	}

}

void HalfFeatures(std::vector<Keypoint>& features)
{
	for(size_t i = 0; i < features.size(); i++)
	{
		features[i].dx /= 2;
		features[i].dy /= 2;
		features[i].scale /= 2;
	}
}

void OrientationAssignment(std::vector<Keypoint>& extrema, std::vector<Keypoint>& features, 
	     const std::vector<Block>& gauss_pyr)
{
	int n = extrema.size();
	double* hist;

	for (int i = 0; i < n; ++i)
	{
		hist = CalculateOrientationHistogram(gauss_pyr[extrema[i].octave*(INTERVALS+3)+extrema[i].interval],
			        extrema[i].x, extrema[i].y, ORI_HIST_BINS, int(ORI_WINDOW_RADIUS*extrema[i].octave_scale+0.5),
					ORI_SIGMA_TIMES*extrema[i].octave_scale);
		for (int j = 0; j < ORI_SMOOTH_TIMES; ++j)
		{
			GaussSmoothOriHist(hist, ORI_HIST_BINS);
		}
		double highest_peak = DominantDirection(hist, ORI_HIST_BINS);

		CalcOriFeatures(extrema[i], features, hist, ORI_HIST_BINS, highest_peak*ORI_PEAK_RATIO);

		delete []hist;
		hist = NULL;
	}

}

double* CalculateOrientationHistogram(const Block& gauss, int x, int y, int bins, int radius, double sigma)
{
	double* hist = new double[bins];

	for (int i = 0; i < bins; ++i)
	{
		hist[i] = 0;
	}

	double mag, ori;
	double weight;

	int bin;
	const double PI2 = 2.0*PI;
	double econs = -1.0/(2.0*sigma*sigma);

	for (int i = -radius; i <= radius; ++i)
	{
		for (int j = -radius; j <= radius; ++j)
		{
			if (CalcGradMagOri(gauss, x+i, y+j, mag, ori))
			{
				weight = exp((i*i+j*j)*econs);

				bin = int(bins*(PI-ori)/PI2+0.5);
				bin = bin<bins ? bin : 0;

				hist[bin] += mag*weight;
			}
		}
	}

	return hist;
}

bool CalcGradMagOri(const Block& gauss, int x, int y, double& mag, double& ori)
{
	if (x > 0 && x < gauss.nXSize-1 && y > 0 && y < gauss.nYSize-1)
	{
		pixel_t* data = gauss.pBuffer;
		int step = gauss.nXSize;

		double dx = DAt(x+1, y)-DAt(x-1, y);
		double dy = DAt(x, y+1)-DAt(x, y-1);

		mag = sqrt(dx*dx + dy*dy);

		ori = atan2(dy, dx);
		return true;
	}
	return false;
}

void GaussSmoothOriHist(double *hist, int n)
{
	double prev = hist[n-1];
	double temp;
	double h0 = hist[0];

	for(int i = 0; i < n; ++i)
	{
		temp = hist[i];
		hist[i] = 0.25*prev+0.5*hist[i]+0.25*(i+1 >= n ? h0 : hist[i+1]);
		prev = temp;
	}

}

double DominantDirection(double *hist, int n)
{
	double maxd = hist[0];

	for(int i = 1; i < n; ++i)
	{
		if(hist[i] > maxd)
		{
			maxd = hist[i];
		}
	}

	return maxd;
}

void CalcOriFeatures(const Keypoint& keypoint, std::vector<Keypoint>& features, 
	     const double *hist, int n, double mag_thr)
{
	double bin, PI2 = PI*2;
	int l, r;

	for (int i = 0; i < n; ++i)
	{
		l = (i == 0) ? n-1 : i-1;
		r = (i+1)%n;

		if (hist[i] > hist[l] && hist[i] > hist[r] && hist[i] >= mag_thr)
		{
			bin = i+Parabola_Interpolate(hist[l], hist[i], hist[r]);
			bin = (bin < 0) ? (bin+n) :(bin >= n ? (bin-n) : bin);
			Keypoint new_key;
			CopyKeypoint(keypoint, new_key);
			new_key.ori = PI2*bin/n - PI;
			features.push_back(new_key);
		}
	}

}

void CopyKeypoint(const Keypoint& src, Keypoint& dst)
{
	dst.dx = src.dx;
	dst.dy = src.dy;

	dst.interval = src.interval;
	dst.octave = src.octave;
	dst.octave_scale = src.octave_scale;
	dst.offset_interval = src.offset_interval;

	dst.offset_x = src.offset_x;
	dst.offset_y = src.offset_y;

	dst.ori = src.ori;
	dst.scale = src.scale;
	dst.val = src.val;
	dst.x = src.x;
	dst.y = src.y;
}

void DescriptorRepresentation(std::vector<Keypoint>& features, 
	      const std::vector<Block>& gauss_pyr, int bins, int width)
{
	double*** hist;

	for (size_t i = 0; i < features.size(); ++i)
	{
		hist = CalculateDescrHist(gauss_pyr[features[i].octave*(INTERVALS+3)+features[i].interval],
			features[i].x, features[i].y, features[i].octave_scale, features[i].ori, bins, width);

		HistToDescriptor(hist, width, bins, features[i]);

		for (int j = 0; j < width; ++j)
		{
			for (int k = 0; k < width; ++k)
			{
				delete [](hist[j][k]);
			}
			delete [](hist[j]);
		}

		delete []hist;
	}

}

double*** CalculateDescrHist(const Block& gauss, int x, int y, 
	               double octave_scale, double ori, int bins, int width)
{
	double*** hist = new double**[width];

	for (int i = 0; i < width; ++i)
	{
		hist[i] = new double*[width];
		for (int j = 0; j < width; ++j)
		{
			hist[i][j] = new double[bins];
		}
	}

	for (int r = 0; r < width; ++r)
	{
		for (int c = 0; c < width; ++c)
		{
			for (int o = 0; o < bins; ++o)
			{
				hist[r][c][o] = 0;
			}
		}
	}

	double cos_ori = cos(ori);
	double sin_ori = sin(ori);

	double sigma = 0.5*width;
	double conste = -1.0/(2*sigma*sigma);

	double PI2 = PI*2;

	double sub_hist_width = DESCR_SCALE_ADJUST*octave_scale;

	int radius = int((sub_hist_width*sqrt(2.0)*(width+1))/2.0+0.5);

	double grad_ori, grad_mag;

	for (int i = -radius; i <= radius; ++i)
	{
		for (int j = -radius; j <= radius; ++j)
		{
			double rot_x = (cos_ori * j - sin_ori * i) / sub_hist_width;
			double rot_y = (sin_ori * j + cos_ori * i) / sub_hist_width;

			double xbin = rot_x+width/2-0.5;
			double ybin = rot_y+width/2-0.5;

			if (xbin > -1.0 && xbin < width && ybin > -1.0 && ybin < width)
			{
				if (CalcGradMagOri(gauss, x+j, y+i, grad_mag, grad_ori))
				{
					grad_ori = PI-grad_ori-ori;
					while(grad_ori < 0.0)
					{
						grad_ori += PI2;
					}
					while(grad_ori >= PI2)
					{
						grad_ori -= PI2;
					}

					double obin = grad_ori*(bins/PI2);

					double weight = exp(conste*(rot_x*rot_x+rot_y*rot_y));

					InterpHistEntry(hist, xbin, ybin, obin, grad_mag*weight, bins, width);
				}
			}
		}
	}

	return hist;
}


void InterpHistEntry(double ***hist, double xbin, double ybin, double obin, double mag, int bins, int d)
{
	double d_r, d_c, d_o, v_r, v_c, v_o;
	double** row, *h;
	int r0, c0, o0, rb, cb, ob, r, c, o;
	r0 = int(floor(ybin));
	c0 = int(floor(xbin));
	o0 = int(floor(obin));
	d_r = ybin-r0;
	d_c = xbin-c0;
	d_o = obin-o0;

	for (r = 0; r <= 1; ++r)
	{
		rb = r0+r;
		if (rb >= 0 && rb < d)
		{
			v_r = mag*((r == 0) ? 1.0-d_r : d_r);
			row = hist[rb];
			for (c = 0; c <= 1; ++c)
			{
				cb = c0+c;
				if (cb >= 0 && cb < d)
				{
					v_c = v_r*((c == 0) ? 1.0-d_c : d_c);
					h = row[cb];
					for (o = 0; o <= 1; ++o)
					{
						ob = (o0+o)%bins;
						v_o = v_c*((o == 0) ? 1.0-d_o : d_o);
						h[ob] += v_o;
					}
				}
			}
		}
	}

}

void HistToDescriptor(double ***hist, int width, int bins, Keypoint& feature)
{
	int int_val, i, r, c, o, k = 0;

	for (r = 0; r < width; ++r)
	{
		for (c = 0; c < width; ++c)
		{
			for (o = 0; o < bins; ++o)
			{
				feature.descriptor[k++] = hist[r][c][o];
			}
		}
	}

	feature.descr_length = k;
	NormalizeDescr(feature);

	for (i = 0; i < k; ++i)
	{
		if (feature.descriptor[i] > DESCR_MAG_THR)
		{
			feature.descriptor[i] = DESCR_MAG_THR;
		}
	}

	NormalizeDescr(feature);

	for (i = 0; i < k; ++i)
	{
		int_val = int(INT_DESCR_FCTR*feature.descriptor[i]);
		feature.descriptor[i] = min(255, int_val);
	}

}

void NormalizeDescr(Keypoint& feat)
{
	double cur, len_inv, len_sq = 0.0;
	int i, d = feat.descr_length;

	for(i = 0; i < d; ++i)
	{
		cur = feat.descriptor[i];
		len_sq += cur*cur;
	}

	len_inv = 1.0 / sqrt(len_sq);

	for(i = 0; i < d; ++i)
	{
		feat.descriptor[i] *= len_inv;
	}

}

bool FeatureCmp(Keypoint& f1, Keypoint& f2)
{
	return f1.scale < f2.scale;
}

kd_node* kdtree_build(Keypoint* features, int n)
{
	kd_node* kd_root;

	if(!features  ||  n <= 0)
	{
		return NULL;
	}

	kd_root = kd_node_init(features, n);
	expand_kd_node_subtree(kd_root);

	return kd_root;
}

kd_node* kd_node_init(Keypoint* features, int n)
{
	kd_node* node;

	node = (kd_node*)malloc(sizeof(kd_node));
	memset(node, 0, sizeof(kd_node));
	node->ki = -1;
	node->features = features;
	node->n = n;

	return node;
}

void kdtree_release(kd_node* kd_root)
{
	if(!kd_root)
	{
		return;
	}

	kdtree_release(kd_root->kd_left);
	kdtree_release(kd_root->kd_right);
	free(kd_root);
}

void expand_kd_node_subtree(kd_node* node)
{
	if(node->n == 1  ||  node->n == 0)
	{
		node->leaf = 1;
		return;
	}

	assign_part_key(node);
	partition_features(node);

	if(node->kd_left)
	{
		expand_kd_node_subtree(node->kd_left);
	}
	if(node->kd_right)
	{
		expand_kd_node_subtree(node->kd_right);
	}
}

void assign_part_key(kd_node* node)
{
	Keypoint* features;
	double kv, x, mean, var, var_max = 0;
	double* tmp;
	int d, n, i, j, ki = 0;

	features = node->features;
	n = node->n;
	d = features[0].descr_length;

	//取最大方差的那一维作为关键的索引
	for(j = 0; j < d; j++)
	{
		mean = var = 0;

		for(i = 0; i < n; i++)
		{
			mean += features[i].descriptor[j];
		}

		mean /= n;

		for(i = 0; i < n; i++)
		{
			x = features[i].descriptor[j] - mean;
			var += x * x;
		}

		var /= n;

		if(var > var_max)
		{
			ki = j;
			var_max = var;
		}
	}

	
	tmp = (double*)calloc(n, sizeof(double));

	for(i = 0; i < n; i++)
	{
		tmp[i] = features[i].descriptor[ki];
	}

	kv = median_select(tmp, n);
	free(tmp);

	node->ki = ki;
	node->kv = kv;
}

double median_select(double* array, int n)
{
	return rank_select(array, n, (n-1)/2);
}

double rank_select(double* array, int n, int r)
{
	double* tmp, med;
	int gr_5, gr_tot, rem_elts, i, j;

	if(n == 1)
	{
		return array[0];
	}

	gr_5 = n / 5;
	gr_tot = (int)ceil(n/5.0);
	rem_elts = n % 5;
	tmp = array;

	for(i = 0; i < gr_5; i++)
	{
		insertion_sort(tmp, 5);
		tmp += 5;
	}

	insertion_sort(tmp, rem_elts);


	tmp = (double*)calloc(gr_tot, sizeof(double));

	for(i = 0, j = 2; i < gr_5; i++, j += 5)
	{
		tmp[i] = array[j];
	}

	if(rem_elts)
	{
		tmp[i++] = array[n - 1 - rem_elts/2];
	}

	med = rank_select(tmp, i, (i-1)/2);
	free(tmp);


	j = partition_array(array, n, med);

	if(r == j)
	{
		return med;
	}
	else if(r<j)
	{
		return rank_select(array, j, r);
	}
	else
	{
		array += j+1;
		return rank_select(array, (n-j-1), (r-j-1));
	}

}

void partition_features(kd_node* node)
{
	Keypoint* features, tmp;
	double kv;
	int n, ki, p, i, j = -1;

	features = node->features;
	n = node->n;
	ki = node->ki;
	kv = node->kv;
	for(i = 0; i < n; i++)
	{
		if(features[i].descriptor[ki] <= kv)
		{
			tmp = features[++j];
			features[j] = features[i];
			features[i] = tmp;
			if(features[j].descriptor[ki] == kv)
				p = j;
		}
	}
	tmp = features[p];
	features[p] = features[j];
	features[j] = tmp;

	if(j == n-1)
	{
		node->leaf = 1;
		return;
	}

	node->kd_left = kd_node_init(features, j+1);
	node->kd_right = kd_node_init(features+(j+1), (n-j-1));
}

void insertion_sort(double* array, int n)
{
	double k;
	int i, j;

	for(i = 1; i < n; i++)
	{
		k = array[i];
		j = i-1;
		while(j >= 0  &&  array[j] > k)
		{
			array[j+1] = array[j];
			j -= 1;
		}
		array[j+1] = k;
	}

}

int partition_array(double* array, int n, double pivot)
{
	double tmp;
	int p, i, j;

	i = -1;

	for(j = 0; j < n; j++)
	{
		if(array[j] <= pivot)
		{
			tmp = array[++i];
			array[i] = array[j];
			array[j] = tmp;
			if(array[i] == pivot)
				p = i;
		}
	}

	array[p] = array[i];
	array[i] = pivot;

	return i;
}

int kdtree_bbf_knn(kd_node* kd_root, Keypoint* feat, int k, Keypoint*** nbrs, int max_nn_chks)
{
	kd_node* expl;
	min_pq* pq;
	Keypoint* tree_feat, ** _nbrs;
	bbf_data* data;
	int i, t = 0, n = 0;

	if(!nbrs  ||  ! feat  ||  ! kd_root)
	{
		return -1;
	}

	_nbrs = (Keypoint**)calloc(k, sizeof(Keypoint*));
	pq = minpq_init();
	minpq_insert(pq, kd_root, 0);
	while(pq->n > 0  &&  t < max_nn_chks)
	{
		expl = (kd_node*)minpq_extract_min(pq);

		if(!expl)
		{
			goto fail;
		}

		expl = explore_to_leaf(expl, feat, pq);

		if(!expl)
		{
			goto fail;
		}

		for(i = 0; i < expl->n; i++)
		{
			tree_feat = &expl->features[i];
			data = (bbf_data*)malloc(sizeof(bbf_data));
			if(!data)
			{
				goto fail;
			}
			data->old_data = tree_feat->feature_data;
			data->d = descr_dist_sq(feat, tree_feat);
			tree_feat->feature_data = data;
			n += insert_into_nbr_array(tree_feat, _nbrs, n, k);
		}

		t++;
	}

	minpq_release(&pq);

	for(i = 0; i < n; i++)
	{
		data = (bbf_data*)_nbrs[i]->feature_data;
		_nbrs[i]->feature_data = data->old_data;
		free(data);
	}

	*nbrs = _nbrs;
	return n;

fail:
	minpq_release(&pq);

	for(i = 0; i < n; i++)
	{
		data = (bbf_data*)_nbrs[i]->feature_data;
		_nbrs[i]->feature_data = data->old_data;
		free(data);
	}

	free(_nbrs);
	*nbrs = NULL;
	return -1;
}

min_pq* minpq_init()
{
	min_pq* pq;
	pq = (min_pq*)malloc(sizeof(min_pq));
	pq->pq_array = (pq_node*)calloc(MINPQ_INIT_NALLOCD, sizeof(pq_node));
	pq->nallocd = MINPQ_INIT_NALLOCD;
	pq->n = 0;

	return pq;
}

int minpq_insert(min_pq* pq, void* data, int key)
{
	int n = pq->n;

	if(pq->nallocd == n)
	{
		pq->nallocd = array_double((void**)&pq->pq_array,
			pq->nallocd,
			sizeof(pq_node));
		if(!pq->nallocd)
		{
			return 1;
		}
	}

	pq->pq_array[n].data = data;
	pq->pq_array[n].key = INT_MAX;
	decrease_pq_node_key(pq->pq_array, pq->n, key);
	pq->n++;

	return 0;
}

void* minpq_extract_min(min_pq* pq)
{
	void* data;

	if(pq->n < 1)
	{
		return NULL;
	}

	data = pq->pq_array[0].data;
	pq->n--;
	pq->pq_array[0] = pq->pq_array[pq->n];
	restore_minpq_order(pq->pq_array, 0, pq->n);

	return data;
}

void minpq_release(min_pq** pq)
{
	if(!pq)
	{
		return;
	}
	if(*pq  &&  (*pq)->pq_array)
	{
		free((*pq)->pq_array);
		free(*pq);
		*pq = NULL;
	}
}

kd_node* explore_to_leaf(kd_node* node, Keypoint* feat, min_pq* min_pq)
{
	kd_node* unexpl, *expl = node;
	double kv;
	int ki;

	while(expl  &&  ! expl->leaf)
	{
		ki = expl->ki;
		kv = expl->kv;

		if(ki >= feat->descr_length)
		{
			return NULL;
		}
		if(feat->descriptor[ki] <= kv)
		{
			unexpl = expl->kd_right;
			expl = expl->kd_left;
		}
		else
		{
			unexpl = expl->kd_left;
			expl = expl->kd_right;
		}

		if(minpq_insert(min_pq, unexpl, (int)abs(kv - feat->descriptor[ki])))
		{
			return NULL;
		}
	}

	return expl;
}

double descr_dist_sq(Keypoint* f1, Keypoint* f2)
{
	double diff, dsq = 0;
	double* descr1, * descr2;
	int i, d;

	d = f1->descr_length;
	if(f2->descr_length != d)
	{
		return DBL_MAX;
	}
	descr1 = f1->descriptor;
	descr2 = f2->descriptor;

	for(i = 0; i < d; i++)
	{
		diff = descr1[i] - descr2[i];
		dsq += diff*diff;
	}
	return dsq;
}

int insert_into_nbr_array(Keypoint* feat, Keypoint** nbrs, int n, int k)
{
	bbf_data*fdata, *ndata;
	double dn, df;
	int i, ret = 0;

	if(n == 0)
	{
		nbrs[0] = feat;
		return 1;
	}

	fdata = (bbf_data*)feat->feature_data;
	df = fdata->d;
	ndata = (bbf_data*)nbrs[n-1]->feature_data;
	dn = ndata->d;
	if(df >= dn)
	{
		if(n == k)
		{
			feat->feature_data = fdata->old_data;
			free(fdata);
			return 0;
		}
		nbrs[n] = feat;
		return 1;
	}


	if(n < k)
	{
		nbrs[n] = nbrs[n-1];
		ret = 1;
	}
	else
	{
		nbrs[n-1]->feature_data = ndata->old_data;
		free(ndata);
	}
	i = n-2;

	while(i >= 0)
	{
		ndata = (bbf_data*)nbrs[i]->feature_data;
		dn = ndata->d;
		if(dn <= df)
		{
			break;
		}
		nbrs[i+1] = nbrs[i];
		i--;
	}

	i++;
	nbrs[i] = feat;

	return ret;
}

int array_double(void** array, int n, int size)
{
	void* tmp;

	tmp = realloc(*array, 2 * n * size);
	if(!tmp)
	{
		if(*array)
		{
			free(*array);
		}
		*array = NULL;
		return 0;
	}
	*array = tmp;
	return n*2;
}

void decrease_pq_node_key(pq_node* pq_array, int i, int key)
{
	pq_node tmp;

	if(key > pq_array[i].key)
	{
		return;
	}

	pq_array[i].key = key;
	while(i > 0  &&  pq_array[i].key < pq_array[parent(i)].key)
	{
		tmp = pq_array[parent(i)];
		pq_array[parent(i)] = pq_array[i];
		pq_array[i] = tmp;
		i = parent(i);
	}
}

void restore_minpq_order(pq_node* pq_array, int i, int n)
{
	pq_node tmp;
	int l, r, min = i;

	l = left(i);
	r = right(i);
	if(l < n)
	{
		if(pq_array[l].key < pq_array[i].key)
		{
			min = l;
		}
	}
	if(r < n)
	{
		if(pq_array[r].key < pq_array[min].key)
		{
			min = r;
		}
	}

	if(min != i)
	{
		tmp = pq_array[min];
		pq_array[min] = pq_array[i];
		pq_array[i] = tmp;
		restore_minpq_order(pq_array, min, n);
	}
}












