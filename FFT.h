#define	NULL 0



//���f���^���`
struct COMPLEX{
	double	R;
	double	I;
};


//���֐����
enum WINDOW_FUNC{RECTANGLE, HUNNING, HUMMING, BLACKMAN, BLACKMAN_HARRIS};


//�e�e�s�N���X
class FFT{
private:
	//�e�[�u��
	struct COMPLEX*		Wn;												//��]���q�e�[�u��
	double*				Window;											//���֐��e�[�u��
	int*				bitReversal;									//�r�b�g�t���e�[�u��



	//�f�[�^
	struct COMPLEX*		DATA;											//�f�[�^

	int					df;												//���g���Ԋu


	//����ϐ�
	int					n;												//����\
	int					SamplingFreq;									//�T���v�����O���g��(44.1kHz��44100)
	WINDOW_FUNC			WindowFunc;										//���֐�

	bool				IsInitialized;									//�������ς݂��ǂ���

	struct COMPLEX		temp;											//�r���v�Z�Ŏg�p


	//�v�Z�⏕�֐�
	//�e�[�u���v�Z
	void				Cal_Rotation();									//��]���q�v�Z
	int					Cal_Window();									//���֐��v�Z
	void				Cal_BitReversal();								//�r�b�g�t���v�Z

	//�o�^�t���C���Z
	void				Butterfly(struct COMPLEX& a, struct COMPLEX& b, int m);

	//��]���q��Z
	void				MulWn(struct COMPLEX& a, int m);


public:
	//�C���^�[�t�F�C�X
	//�������֌W
						FFT();											//�R���X�g���N�^
						~FFT();											//�f�X�g���N�^

	//�������֐�
	int					Init(int k, WINDOW_FUNC WF = HUNNING, int S_Freq = 44100);

	void				Reset();										//���Z�b�g�֐�

	void				GetBitTable(int* a);							//�r�b�g�t���e�[�u���擾

	//�e�e�s�v�Z�֐�
	void				Calculate(struct COMPLEX& data, double* DFT = NULL);
};


class EFFECTOR : public FFT{
private:
	struct COMPLEX*		S_DATA;											//IFFT��̐M���f�[�^

	//�⏕�v�Z�֐�
	void				MulWn_I(struct COMPLEX& a, int m);				//IFFT���̉�]���q��Z
	//IFFT���̃o�^�t���C���Z
	void				Butterfly_I(struct COMPLEX& a, struct COMPLEX& b, int m);

public:
	//�h�e�e�s�v�Z�֐�
	void				IFFT(struct COMPLEX& s_data);
};


//�~����
const double	PI		= 3.14159265358979;



