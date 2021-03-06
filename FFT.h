#define	NULL 0



//¡f^ðè`
struct COMPLEX{
	double	R;
	double	I;
};


//Öðñ
enum WINDOW_FUNC{RECTANGLE, HUNNING, HUMMING, BLACKMAN, BLACKMAN_HARRIS};


//eesNX
class FFT{
private:
	//e[u
	struct COMPLEX*		Wn;												//ñ]öqe[u
	double*				Window;											//Öe[u
	int*				bitReversal;									//rbgte[u



	//f[^
	struct COMPLEX*		DATA;											//f[^

	int					df;												//ügÔu


	//§äÏ
	int					n;												//ªð\
	int					SamplingFreq;									//TvOüg(44.1kHz¨44100)
	WINDOW_FUNC			WindowFunc;										//Ö

	bool				IsInitialized;									//ú»ÏÝ©Ç¤©

	struct COMPLEX		temp;											//rvZÅgp


	//vZâÖ
	//e[uvZ
	void				Cal_Rotation();									//ñ]öqvZ
	int					Cal_Window();									//ÖvZ
	void				Cal_BitReversal();								//rbgtvZ

	//o^tCZ
	void				Butterfly(struct COMPLEX& a, struct COMPLEX& b, int m);

	//ñ]öqæZ
	void				MulWn(struct COMPLEX& a, int m);


public:
	//C^[tFCX
	//ú»ÖW
						FFT();											//RXgN^
						~FFT();											//fXgN^

	//ú»Ö
	int					Init(int k, WINDOW_FUNC WF = HUNNING, int S_Freq = 44100);

	void				Reset();										//ZbgÖ

	void				GetBitTable(int* a);							//rbgte[uæ¾

	//eesvZÖ
	void				Calculate(struct COMPLEX& data, double* DFT = NULL);
};


class EFFECTOR : public FFT{
private:
	struct COMPLEX*		S_DATA;											//IFFTãÌMf[^

	//âvZÖ
	void				MulWn_I(struct COMPLEX& a, int m);				//IFFTÌñ]öqæZ
	//IFFTÌo^tCZ
	void				Butterfly_I(struct COMPLEX& a, struct COMPLEX& b, int m);

public:
	//heesvZÖ
	void				IFFT(struct COMPLEX& s_data);
};


//~ü¦
const double	PI		= 3.14159265358979;



