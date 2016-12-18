

#include"FFT.h"
#include<math.h>
#include<windows.h>
#include<stdio.h>



//�R���X�g���N�^
FFT::FFT()
{
	Wn				= NULL;				//��]���q�e�[�u���ւ̃|�C���^��������
	Window			= NULL;				//���֐��e�[�u���ւ̃|�C���^��������
	bitReversal		= NULL;				//�r�b�g�t���e�[�u���ւ̃|�C���^��������


	DATA			= NULL;				//�M���f�[�^�ւ̃|�C���^��������


	IsInitialized	= FALSE;			//�����������
}


//�f�X�g���N�^
FFT::~FFT()
{
	//�������ς݂Ȃ烁�����J��
	if(IsInitialized)
	{
		//DFT�l�f�[�^���J��
		delete[]	DATA;
		DATA		= NULL;

		//��]���q�e�[�u�����J��
		delete[]	Wn;
		Wn			= NULL;

		//�r�b�g�t���e�[�u�����J��
		delete[]	bitReversal;
		bitReversal	= NULL;
	}
}


//�������֐�
int		FFT::Init(int k, WINDOW_FUNC WF, int S_Freq)
{
	if(IsInitialized)
		return 0;

	//�e��ݒ�
	n				= k;					//����\�ݒ�
	WindowFunc		= WF;					//���֐��ݒ�
	SamplingFreq	= S_Freq;				//�T���v�����O���g���ݒ�


	//����\��2�ׂ̂���łȂ���΃G���[
	if(k % 2)
		return -2;


	//�������m��

	//�f�[�^�z��
	if((DATA = new struct COMPLEX[n]) == NULL)
		return -1;

	//�f�[�^�e�[�u��

	//��]���q�e�[�u��
	if((Wn = new struct COMPLEX[n]) == NULL)
		return -1;

	//���֐��e�[�u��
	if((Window = new double[n]) == NULL)
		return -1;

	//�r�b�g�t���e�[�u��
	if((bitReversal = new int[n]) == NULL)
		return -1;

	//���֐��e�[�u���v�Z
	if(Cal_Window())
		return -3;

	//��]���q�e�[�u���v�Z
	Cal_Rotation();

	//�r�b�g�t���v�Z
	Cal_BitReversal();


	IsInitialized	= TRUE;					//�������ς�

	//���g���Ԋu�v�Z
	df				= SamplingFreq / n;

	return 0;
}


int FFT::Cal_Window()
{
	switch(WindowFunc)
	{
	//��`��
	case RECTANGLE:
		for(int i = 0; i < n; i++)
			Window[i] = 1;
		break;

	//�n�j���O��
	case HUNNING:
		for(int i = 0; i < n; i++)
			Window[i] = 0.5 - 0.5 * cos((2 * PI * (double)i) / (double)n);
		break;

	//�n�~���O��
	case HUMMING:
		for(int i = 0; i < n; i++)
			Window[i] = 0.54 - 0.46 * cos((2 * PI * (double)i) / (double)n);
		break;

	//�u���b�N�}����
	case BLACKMAN:
		for(int i = 0; i < n; i++)
			Window[i] = 0.42 - 0.5 * cos((2 * PI * (double)i) / (double)n) + 0.08 * cos((4 * PI * (double)i) / (double)n);
		break;

	//�u���b�N�}���n���X��
	case BLACKMAN_HARRIS:
		for(int i = 0; i < n; i++)
			Window[i] = 0.35875 - 0.48829 * cos((2 * PI * (double)i) / (double)n) + 0.14128 * cos((4 * PI * (double)i) / (double)n) - 0.01168 * cos((6 * PI * (double)i) / (double)n);
		break;

	default:
		return 1;
	}
	return 0;
}


//��]���q�e�[�u������
void	FFT::Cal_Rotation()
{
	for(int i = 0; i < (n / 2); i++)
	{
		Wn[i].R	= cos((2 * PI * (double)i) / (double)n);
		Wn[i].I	= -sin((2 * PI * (double)i) / (double)n);
	}
}

//�r�b�g�t���v�Z
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



//��]���q��Z
void	FFT::MulWn(struct COMPLEX &a, int m)
{
	temp.R	= a.R * Wn[m].R - a.I * Wn[m].I;
	a.I		= a.R * Wn[m].I + a.I * Wn[m].R;

	a.R		= temp.R;
}


//�o�^�t���C���Z
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
	int			i, j, k, l;				//���[�v����
	int			NOD;					//������


	//�f�[�^���R�s�[
	for(i = 0; i < n; i++)
	{
		DATA[i].R	= (&data)[i].R;
		DATA[i].I	= (&data)[i].I;
	}

	//���|��
	for(i = 0; i < n; i++)
	{
		DATA[i].R	*= Window[i];
		DATA[i].I	*= Window[i];
	}

	//FFT�v�Z�i���g���Ԉ����j
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


//�O������r�b�g�t���e�[�u�����擾����֐�
void	FFT::GetBitTable(int* a)
{
	for(int i = 0; i < n; i++)
		a[i]	= bitReversal[i];
}



void	FFT::Reset()
{
		//�������ς݂Ȃ烁�����J��
	if(IsInitialized)
	{
		//�f�[�^���J��
		delete[] DATA;
		DATA			= NULL;

		//��]���q�e�[�u�����J��
		delete[] Wn;
		Wn				= NULL;

		//�r�b�g�t���e�[�u�����J��
		delete[] bitReversal;
		bitReversal		= NULL;


		IsInitialized	= FALSE;
	}
}

