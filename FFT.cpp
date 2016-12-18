

#include"FFT.h"
#include<math.h>
#include<windows.h>
#include<stdio.h>



//コンストラクタ
FFT::FFT()
{
	Wn				= NULL;				//回転因子テーブルへのポインタを初期化
	Window			= NULL;				//窓関数テーブルへのポインタを初期化
	bitReversal		= NULL;				//ビット逆順テーブルへのポインタを初期化


	DATA			= NULL;				//信号データへのポインタを初期化


	IsInitialized	= FALSE;			//未初期化状態
}


//デストラクタ
FFT::~FFT()
{
	//初期化済みならメモリ開放
	if(IsInitialized)
	{
		//DFT値データを開放
		delete[]	DATA;
		DATA		= NULL;

		//回転因子テーブルを開放
		delete[]	Wn;
		Wn			= NULL;

		//ビット逆順テーブルを開放
		delete[]	bitReversal;
		bitReversal	= NULL;
	}
}


//初期化関数
int		FFT::Init(int k, WINDOW_FUNC WF, int S_Freq)
{
	if(IsInitialized)
		return 0;

	//各種設定
	n				= k;					//分解能設定
	WindowFunc		= WF;					//窓関数設定
	SamplingFreq	= S_Freq;				//サンプリング周波数設定


	//分解能が2のべき乗でなければエラー
	if(k % 2)
		return -2;


	//メモリ確保

	//データ配列
	if((DATA = new struct COMPLEX[n]) == NULL)
		return -1;

	//データテーブル

	//回転因子テーブル
	if((Wn = new struct COMPLEX[n]) == NULL)
		return -1;

	//窓関数テーブル
	if((Window = new double[n]) == NULL)
		return -1;

	//ビット逆順テーブル
	if((bitReversal = new int[n]) == NULL)
		return -1;

	//窓関数テーブル計算
	if(Cal_Window())
		return -3;

	//回転因子テーブル計算
	Cal_Rotation();

	//ビット逆順計算
	Cal_BitReversal();


	IsInitialized	= TRUE;					//初期化済み

	//周波数間隔計算
	df				= SamplingFreq / n;

	return 0;
}


int FFT::Cal_Window()
{
	switch(WindowFunc)
	{
	//矩形窓
	case RECTANGLE:
		for(int i = 0; i < n; i++)
			Window[i] = 1;
		break;

	//ハニング窓
	case HUNNING:
		for(int i = 0; i < n; i++)
			Window[i] = 0.5 - 0.5 * cos((2 * PI * (double)i) / (double)n);
		break;

	//ハミング窓
	case HUMMING:
		for(int i = 0; i < n; i++)
			Window[i] = 0.54 - 0.46 * cos((2 * PI * (double)i) / (double)n);
		break;

	//ブラックマン窓
	case BLACKMAN:
		for(int i = 0; i < n; i++)
			Window[i] = 0.42 - 0.5 * cos((2 * PI * (double)i) / (double)n) + 0.08 * cos((4 * PI * (double)i) / (double)n);
		break;

	//ブラックマンハリス窓
	case BLACKMAN_HARRIS:
		for(int i = 0; i < n; i++)
			Window[i] = 0.35875 - 0.48829 * cos((2 * PI * (double)i) / (double)n) + 0.14128 * cos((4 * PI * (double)i) / (double)n) - 0.01168 * cos((6 * PI * (double)i) / (double)n);
		break;

	default:
		return 1;
	}
	return 0;
}


//回転因子テーブル生成
void	FFT::Cal_Rotation()
{
	for(int i = 0; i < (n / 2); i++)
	{
		Wn[i].R	= cos((2 * PI * (double)i) / (double)n);
		Wn[i].I	= -sin((2 * PI * (double)i) / (double)n);
	}
}

//ビット逆順計算
void	FFT::Cal_BitReversal()
{
	int		a		= 0;

	bitReversal[0]	= 0;

	for(int i = 1; i < n; i++)
	{
		a	^= (n / 2);
		for(int j = 2; j <= n; j *= 2)
		{
			if(!(i % j))
				a	^= (n / (2 * j));
		}
		bitReversal[i]	= a;
	}
}



//回転因子乗算
void	FFT::MulWn(struct COMPLEX &a, int m)
{
	temp.R	= a.R * Wn[m].R - a.I * Wn[m].I;
	a.I		= a.R * Wn[m].I + a.I * Wn[m].R;

	a.R		= temp.R;
}


//バタフライ演算
void	FFT::Butterfly(struct COMPLEX& a, struct COMPLEX& b, int m)
{
	temp.R		= a.R + b.R;
	temp.I		= a.I + b.I;

	b.R			= a.R - b.R;
	b.I			= a.I - b.I;

	a.R			= temp.R;
	a.I			= temp.I;

	MulWn(b, m);
}


void	FFT::Calculate(struct COMPLEX& data, double* DFT)
{
	int			i, j, k, l;				//ループ制御
	int			NOD;					//分割数


	//データをコピー
	for(i = 0; i < n; i++)
	{
		DATA[i].R	= (&data)[i].R;
		DATA[i].I	= (&data)[i].I;
	}

	//窓掛け
	for(i = 0; i < n; i++)
	{
		DATA[i].R	*= Window[i];
		DATA[i].I	*= Window[i];
	}

	//FFT計算（周波数間引き）
	for(NOD = 1; NOD < n; NOD *= 2)
	{
		l = 0;
		for(j = 0; j < NOD; j ++)
		{
			for(k = 0; k < (n / (2 * NOD)); k++)
				Butterfly(DATA[k + l], DATA[k + l + (n / (2 * NOD))], k * NOD);

			l += (n / NOD);
		}
	}

	for(int i = 0; i < n; i++)
	{
		DATA[i].R /= (double)n;
		DATA[i].I /= (double)n;
	}

	if(DFT == NULL)
		return ;

	for(int i = 0; i < n; i++)
		DFT[i]	= 2 * sqrt(DATA[bitReversal[i]].R * DATA[bitReversal[i]].R + DATA[bitReversal[i]].I * DATA[bitReversal[i]].I);

	return ;
}


//外部からビット逆順テーブルを取得する関数
void	FFT::GetBitTable(int* a)
{
	for(int i = 0; i < n; i++)
		a[i]	= bitReversal[i];
}



void	FFT::Reset()
{
		//初期化済みならメモリ開放
	if(IsInitialized)
	{
		//データを開放
		delete[] DATA;
		DATA			= NULL;

		//回転因子テーブルを開放
		delete[] Wn;
		Wn				= NULL;

		//ビット逆順テーブルを開放
		delete[] bitReversal;
		bitReversal		= NULL;


		IsInitialized	= FALSE;
	}
}

