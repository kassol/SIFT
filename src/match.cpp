#include "match.h"
#include "matParamEstimator.h"
#include "ransac.h"
#include <algorithm>




match::match(char* l, char* r)
{
	m_szPathNameL = l;
	m_szPathNameR = r;
	CoInitialize(NULL);
	m_pImage = NULL;
	m_pImage2 = NULL;
	::CoCreateInstance(CLSID_ImageDriver, NULL, CLSCTX_ALL, IID_IImage, (void**)&m_pImage);
	::CoCreateInstance(CLSID_ImageDriver, NULL, CLSCTX_ALL, IID_IImage, (void**)&m_pImage2);
}


match::~match(void)
{
	m_pImage->Release();
	::CoUninitialize();
}

void match::domatch(std::vector<SamePoint>& resultData)
{
	const int nBlockSize = 512;
	int scale = 1;
	int nx1, ny1, nx2, ny2, nband1, nband2;
	uchar *pBuf = NULL;
	m_pImage->Open(_bstr_t(m_szPathNameL), modeRead);
	m_pImage->GetCols(&nx1);
	m_pImage->GetRows(&ny1);
	m_pImage->GetBandNum(&nband1);
	m_pImage->Close();
	m_pImage->Open(_bstr_t(m_szPathNameR), modeRead);
	m_pImage->GetCols(&nx2);
	m_pImage->GetRows(&ny2);
	m_pImage->GetBandNum(&nband2);
	m_pImage->Close();
	int mincr = min(min(min(nx1, ny1), nx2),ny2);
	while(mincr/scale >1024)
	{
		scale *= 2;
	}
	int nxsize1 = nx1/scale;
	int nysize1 = ny1/scale;
	int nxsize2 = nx2/scale;
	int nysize2 = ny2/scale;

	m_pImage->Open(_bstr_t(m_szPathNameL), modeRead);
	pBuf = new uchar[nxsize1*nysize1*nband1];
	m_pImage->ReadImg(0, 0, nx1, ny1, pBuf, nxsize1, nysize1, nband1, 0, 0, nxsize1, nysize1, -1, 0);
	m_pImage->Close();
	m_pImage->CreateImg(_bstr_t("templ.tif"), modeCreate, nxsize1, nysize1, Pixel_Byte, nband1, BIL, 0, 0, 1);
	m_pImage->WriteImg(0, 0, nxsize1, nysize1, pBuf, nxsize1, nysize1, nband1, 0, 0, nxsize1, nysize1, -1, 0);
	m_pImage->Close();
	delete [] pBuf;
	

	m_pImage->Open(_bstr_t(m_szPathNameR), modeRead);
	pBuf = new uchar[nxsize2*nysize2*nband2];
	m_pImage->ReadImg(0, 0, nx2, ny2, pBuf, nxsize2, nysize2, nband2, 0, 0, nxsize2, nysize2, -1, 0);
	m_pImage->Close();
	m_pImage->CreateImg(_bstr_t("tempr.tif"), modeCreate, nxsize2, nysize2, Pixel_Byte, nband2, BIL, 0, 0, 1);
	m_pImage->WriteImg(0, 0, nxsize2, nysize2, pBuf, nxsize2, nysize2, nband2, 0, 0, nxsize2, nysize2, -1, 0);
	m_pImage->Close();
	delete []pBuf;
	pBuf = NULL;

	sl.fetchFeatures("templ.tif");
	sr.fetchFeatures("tempr.tif");

	/*sl.fetchFeatures(m_szPathNameL);
	sr.fetchFeatures(m_szPathNameR);*/

	Keypoint** nbrs;
	kd_node* kd_root;
	int nsize = sl.m_listKeyPoint.size();
	int nsize2 = sr.m_listKeyPoint.size();
	Keypoint* feat1 = (Keypoint*)malloc(nsize*sizeof(Keypoint));
	std::list<Keypoint>::iterator temIte = sl.m_listKeyPoint.begin();
	int i = 0;
	while(temIte != sl.m_listKeyPoint.end())
	{
		feat1[i] = *temIte;
		++i;
		++temIte;
	}
	sl.m_listKeyPoint.clear();
	Keypoint* feat2 = (Keypoint*)malloc(nsize2*sizeof(Keypoint));
	temIte = sr.m_listKeyPoint.begin();
	i = 0;
	while(temIte != sr.m_listKeyPoint.end())
	{
		feat2[i] = *temIte;
		++i;
		++temIte;
	}
	sr.m_listKeyPoint.clear();

	kd_root = kdtree_build(feat2, nsize2);
	int k = 0;
	double d0, d1;
	Keypoint* feat;
	int matchnum = 0;

	std::vector<SamePoint> sp;
	for (i = 0; i < nsize; ++i)
	{
		feat = feat1+i;
		k = kdtree_bbf_knn(kd_root, feat, 2, &nbrs, KDTREE_BBF_MAX_NN_CHKS);
		if (k == 2)
		{
			d0 = descr_dist_sq(feat, nbrs[0]);
			d1 = descr_dist_sq(feat, nbrs[1]);
			if (d0 < d1*NN_SQ_DIST_RATIO_THR)
			{
				sp.push_back(SamePoint(feat->dx, feat->dy, nbrs[0]->dx, nbrs[0]->dy));
				++matchnum;
			}
		}
		free(nbrs);
	}
	kdtree_release(kd_root);

	free(feat1);
	free(feat2);
	feat1 = NULL;
	feat2 = NULL;


	std::vector<double> matParameters;
	MatParamEstimator mpEstimator(5);
	int numForEstimate = 3;

	mpEstimator.leastSquaresEstimate(sp, matParameters);

	double usedData = Ransac<SamePoint, double>::compute(matParameters, &mpEstimator, sp, numForEstimate, 0.9, 0.1, resultData);
	sp.swap(std::vector<SamePoint>());
	resultData.swap(std::vector<SamePoint>());
	
	Pt lLeftTop(0, 0);
	Pt lRightBottom(nx1, ny1);
	double a = matParameters[0];
	double b = matParameters[1];
	double c = matParameters[2]*scale;
	double d = matParameters[3];
	double e = matParameters[4];
	double f = matParameters[5]*scale;

	std::list<Block> listDataBlock;
	for (int y = 0; y < ny1;)
	{
		if (y+nBlockSize < ny1)
		{
			for (int x = 0; x < nx1;)
			{
				if (x+nBlockSize < nx1)
				{
					Rec rect(a*x+b*y+c, a*(x+nBlockSize)+b*y+c, d*x+e*y+f, d*x+e*(y+nBlockSize)+f);
					if (rect.Intersects(Rec(0, nx2, 0, ny2)))
					{
						listDataBlock.push_back(Block(NULL, x, y, nBlockSize, nBlockSize));
					}
					x += nBlockSize;
				}
				else
				{
					Rec rect(a*x+b*y+c, a*nx1+b*y+c, d*x+e*y+f, d*x+e*(y+nBlockSize)+f);
					if (rect.Intersects(Rec(0, nx2, 0, ny2)))
					{
						listDataBlock.push_back(Block(NULL, x, y, nx1-x, nBlockSize));
					}
					x = nx1;
				}
			}
			y += nBlockSize;
		}
		else
		{
			for (int x = 0; x < nx1;)
			{
				if (x+nBlockSize < nx1)
				{
					Rec rect(a*x+b*y+c, a*(x+nBlockSize)+b*y+c, d*x+e*y+f, d*x+e*ny1+f);
					if (rect.Intersects(Rec(0, nx2, 0, ny2)))
					{
						listDataBlock.push_back(Block(NULL, x, y, nBlockSize, ny1-y));
					}
					x += nBlockSize;
				}
				else
				{
					Rec rect(a*x+b*y+c, a*nx1+b*y+c, d*x+e*y+f, d*x+e*ny1+f);
					if (rect.Intersects(Rec(0, nx2, 0, ny2)))
					{
						listDataBlock.push_back(Block(NULL, x, y, nx1-x, ny1-y));
					}
					x = nx1;
				}
			}
			y = ny1;
		}
	}
	int nBlockNumx = (nx2+nBlockSize-1)/nBlockSize;
	int nBlockNumy = (ny2+nBlockSize-1)/nBlockSize;
	int nBlockNum = nBlockNumx*nBlockNumy;
	std::vector<RBTree> vecKDTree(nBlockNum);
	std::list<Block>::iterator blockIte = listDataBlock.begin();
	m_pImage->Open(_bstr_t(m_szPathNameL), modeRead);
	m_pImage2->Open(_bstr_t(m_szPathNameR), modeRead);
	while(blockIte != listDataBlock.end())
	{
		pBuf = new uchar[blockIte->nXSize*blockIte->nYSize*nband1];
		m_pImage->ReadImg(blockIte->nXOrigin, blockIte->nYOrigin, 
			blockIte->nXOrigin+blockIte->nXSize, blockIte->nYOrigin+blockIte->nYSize, pBuf,
			blockIte->nXSize, blockIte->nYSize, nband1, 0, 0,
			blockIte->nXSize, blockIte->nYSize, -1, 0);
		pixel_t* p = new pixel_t[blockIte->nXSize*blockIte->nYSize];
		for (int y = 0; y < blockIte->nYSize; ++y)
		{
			for (int x = 0, m = 0; x < blockIte->nXSize*nband1; x += nband1, ++m)
			{
				double sum = 0;
				for (int n = 0; n < nband1; ++n)
				{
					sum += pBuf[y*blockIte->nXSize*nband1+x+n];
				}
				p[y*blockIte->nXSize+m] = sum/(nband1*225.0);
			}
		}
		delete []pBuf;
		pBuf = NULL;
		std::vector<Keypoint> feature;
		sift(p, feature, blockIte->nXSize, blockIte->nYSize);
		p = NULL;
		std::vector<Keypoint>::iterator feaIte = feature.begin();
		while(feaIte != feature.end())
		{
			feaIte->dx += blockIte->nXOrigin;
			feaIte->dy += blockIte->nYOrigin;
			int calx = int(feaIte->dx*a+feaIte->dy*b+c);
			int caly = int(feaIte->dx*d+feaIte->dy*e+f);
			int idx = calx/nBlockSize;
			int idy = caly/nBlockSize;
			int nBlockIndex = idy*nBlockNumx+idx;
			if (vecKDTree[nBlockIndex].num != 0 && vecKDTree[nBlockIndex].num < 50)
			{
				vecKDTree[nBlockIndex].num == 0;
				++feaIte;
				continue;
			}
			if (vecKDTree[nBlockIndex].node == NULL)
			{
				int xo = idx*nBlockSize;
				int yo = idy*nBlockSize;
				int xsize = nBlockSize;
				int ysize = nBlockSize;
				if (idx == nBlockNumx-1)
				{
					xsize = nx2%nBlockSize;
				}
				if (idy == nBlockNumy-1)
				{
					ysize = ny2%nBlockSize;
				}
				pBuf = new uchar[xsize*ysize*nband2];
				m_pImage2->ReadImg(xo, yo, xo+xsize, yo+ysize, pBuf, xsize, ysize, nband2, 0, 0, xsize, ysize, -1, 0);
				p = new pixel_t[xsize*ysize];
				for(int y = 0; y < ysize; ++y)
				{
					for (int x = 0, m = 0; x < xsize*nband2; x += nband2, ++m)
					{
						double sum = 0;
						for (int n = 0; n < nband2; ++n)
						{
							sum += pBuf[y*xsize*nband2+x+n];
						}
						p[y*xsize+m] = sum/(nband2*255.0);
					}
				}
				delete []pBuf;
				pBuf = NULL;
				std::vector<Keypoint> feature2;
				sift(p, feature2, xsize, ysize);
				p = NULL;
				int nf2 = feature2.size();
				vecKDTree[nBlockIndex].num = nf2;
				if (nf2 < 50)
				{
					++feaIte;
					continue;
				}
				feat2 = (Keypoint*)malloc(nf2*sizeof(Keypoint));
				std::vector<Keypoint>::iterator kIte2 = feature2.begin();
				i = 0;
				while(kIte2 != feature2.end())
				{
					kIte2->dx += xo;
					kIte2->dy += yo;
					feat2[i] = *kIte2;
					++i;
					++kIte2;
				}
				feature2.swap(std::vector<Keypoint>());
				kd_root = kdtree_build(feat2, nf2);
				vecKDTree[nBlockIndex].node = kd_root;
				vecKDTree[nBlockIndex].feature = feat2;
			}
			k = kdtree_bbf_knn(vecKDTree[nBlockIndex].node, &(*feaIte), 2, &nbrs, KDTREE_BBF_MAX_NN_CHKS);
			if (k == 2)
			{
				d0 = descr_dist_sq(&(*feaIte), nbrs[0]);
				d1 = descr_dist_sq(&(*feaIte), nbrs[1]);
				if (d0 < d1*NN_SQ_DIST_RATIO_THR)
				{
					sp.push_back(SamePoint(feaIte->dx, feaIte->dy, nbrs[0]->dx, nbrs[0]->dy));
				}
			}
			free(nbrs);
			++feaIte;
		}
		feature.swap(std::vector<Keypoint>());
		std::vector<RBTree>::iterator kdIte = vecKDTree.begin();
		while(kdIte != vecKDTree.end())
		{
			if (kdIte->node != NULL)
			{
				free(kdIte->node);
				kdIte->node = NULL;
				free(kdIte->feature);
				kdIte->feature = NULL;
			}
			++kdIte;
		}
		++blockIte;
	}
	m_pImage->Close();
	m_pImage2->Close();

	mpEstimator.leastSquaresEstimate(sp, matParameters);

	usedData = Ransac<SamePoint, double>::compute(matParameters, &mpEstimator, sp, numForEstimate, 0.9, 0.1, resultData);
	sp.swap(std::vector<SamePoint>());
	resultData.swap(std::vector<SamePoint>());

	a = matParameters[0];
	b = matParameters[1];
	c = matParameters[2];
	d = matParameters[3];
	e = matParameters[4];
	f = matParameters[5];


	//mosaic
	Pt calrLeftTop;
	Pt calrRightBottom;
	calrLeftTop.x = lLeftTop.x*a+lLeftTop.y*b+c;
	calrLeftTop.y = lLeftTop.x*d+lLeftTop.y*e+f;
	calrRightBottom.x = lRightBottom.x*a+lRightBottom.y*b+c;
	calrRightBottom.y = lRightBottom.x*d+lRightBottom.y*e+f;
	Rec calRight(calrLeftTop.x, calrRightBottom.x, calrLeftTop.y, calrRightBottom.y);
	
	Rec Right(0, nx2, 0, ny2);
	Rec resultRectR = calRight.Intersected(Right);
	resultRectR.extend(nx2, ny2);
	Pt callLeftTop;
	Pt callRightBottom;
	callLeftTop.y = (d/a*resultRectR.left-resultRectR.top+f-c*d/a)/(b*d/a-e);
	callLeftTop.x = (resultRectR.left-b*callLeftTop.y-c)/a;
	callRightBottom.y = (d/a*resultRectR.right-resultRectR.bottom+f-c*d/a)/(b*d/a-e);
	callRightBottom.x = (resultRectR.right-b*callRightBottom.y-c)/a;
	Rec resultRectL(callLeftTop.x, callRightBottom.x, callLeftTop.y, callRightBottom.y);
	resultRectL.extend(nx1, ny1);


	int newXsize, newYsize;
	newXsize = (int)max(nx2, calrRightBottom.x);
	newYsize = (int)max(ny2, calrRightBottom.y);
	
	m_pImage->CreateImg(_bstr_t("result.tif"), modeCreate, newXsize, newYsize, Pixel_Byte, nband1, BIL, 0, 0, 1);
	IImage* pImage = NULL;
	::CoCreateInstance(CLSID_ImageDriver, NULL, CLSCTX_ALL, IID_IImage, (void**)&pImage);
	pImage->Open(_bstr_t(m_szPathNameR), modeRead);
	pBuf = new uchar[nx2*nBlockSize*nband2];
	for (int i = 0; i < ny2;)
	{
		if (i + nBlockSize < ny2)
		{
			pImage->ReadImg(0, i, nx2, i+nBlockSize, pBuf, nx2, nBlockSize, nband2, 0, 0, nx2, nBlockSize, -1, 0);
			m_pImage->WriteImg(0, i, nx2, i+nBlockSize, pBuf, nx2, nBlockSize, nband2, 0, 0, nx2, nBlockSize, -1, 0);
			i += nBlockSize;
		}
		else
		{
			pImage->ReadImg(0, i, nx2, ny2, pBuf, nx2, nBlockSize, nband2, 0, 0, nx2, ny2-i, -1, 0);
			m_pImage->WriteImg(0, i, nx2, ny2, pBuf, nx2, nBlockSize, nband2, 0, 0, nx2, ny2-i, -1, 0);
			i = ny2;
		}
	}
	pImage->Close();
	delete []pBuf;
	pImage->Open(_bstr_t(m_szPathNameL), modeRead);
	pBuf = new uchar[nx1*nBlockSize*nband1];
	for (int i = 0; i < ny1;)
	{
		if (i+nBlockSize < ny1)
		{
			pImage->ReadImg(0, i, nx1, i+nBlockSize, pBuf, nx1, nBlockSize, nband1, 0, 0, nx1, nBlockSize, -1, 0);
			m_pImage->WriteImg((int)calrLeftTop.x, (int)calrLeftTop.y+i, (int)calrRightBottom.x, (int)calrLeftTop.y+i+nBlockSize, pBuf, nx1, nBlockSize, nband1, 0, 0, nx1, nBlockSize, -1, 0);
			i += nBlockSize;
		}
		else
		{
			pImage->ReadImg(0, i, nx1, ny1, pBuf, nx1, nBlockSize, nband1, 0, 0, nx1, ny1-i, -1, 0);
			m_pImage->WriteImg((int)calrLeftTop.x, (int)calrLeftTop.y+i, (int)calrRightBottom.x, (int)calrLeftTop.y+ny1, pBuf, nx1, nBlockSize, nband1, 0, 0, nx1, ny1-i, -1, 0);
			i = ny1;
		}
	}
	pImage->Close();
	delete []pBuf;
	pBuf = NULL;
	m_pImage->Close();
}
