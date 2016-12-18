#define	NULL 0



//複素数型を定義
struct COMPLEX{
	double	R;
	double	I;
};


//窓関数を列挙
enum WINDOW_FUNC{RECTANGLE, HUNNING, HUMMING, BLACKMAN, BLACKMAN_HARRIS};


//ＦＦＴクラス
class FFT{
private:
	//テーブル
	struct COMPLEX*		Wn;												//回転因子テーブル
	double*				Window;											//窓関数テーブル
	int*				bitReversal;									//ビット逆順テーブル



	//データ
	struct COMPLEX*		DATA;											//データ

	int					df;												//周波数間隔


	//制御変数
	int					n;												//分解能
	int					SamplingFreq;									//サンプリング周波数(44.1kHz→44100)
	WINDOW_FUNC			WindowFunc;										//窓関数

	bool				IsInitialized;									//初期化済みかどうか

	struct COMPLEX		temp;											//途中計算で使用


	//計算補助関数
	//テーブル計算
	void				Cal_Rotation();									//回転因子計算
	int					Cal_Window();									//窓関数計算
	void				Cal_BitReversal();								//ビット逆順計算

	//バタフライ演算
	void				Butterfly(struct COMPLEX& a, struct COMPLEX& b, int m);

	//回転因子乗算
	void				MulWn(struct COMPLEX& a, int m);


public:
	//インターフェイス
	//初期化関係
						FFT();											//コンストラクタ
						~FFT();											//デストラクタ

	//初期化関数
	int					Init(int k, WINDOW_FUNC WF = HUNNING, int S_Freq = 44100);

	void				Reset();										//リセット関数

	void				GetBitTable(int* a);							//ビット逆順テーブル取得

	//ＦＦＴ計算関数
	void				Calculate(struct COMPLEX& data, double* DFT = NULL);
};


class EFFECTOR : public FFT{
private:
	struct COMPLEX*		S_DATA;											//IFFT後の信号データ

	//補助計算関数
	void				MulWn_I(struct COMPLEX& a, int m);				//IFFT時の回転因子乗算
	//IFFT時のバタフライ演算
	void				Butterfly_I(struct COMPLEX& a, struct COMPLEX& b, int m);

public:
	//ＩＦＦＴ計算関数
	void				IFFT(struct COMPLEX& s_data);
};


//円周率
const double	PI		= 3.14159265358979;



