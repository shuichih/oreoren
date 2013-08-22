#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bmpexporter.h"

#pragma pack(push,2)

typedef struct tagBITMAPFILEHEADER {
  unsigned short bfType;
  unsigned long  bfSize;
  unsigned short bfReserved1;
  unsigned short bfReserved2;
  unsigned long  bfOffBits;
} BITMAPFILEHEADER;

typedef struct tagBITMAPINFOHEADER{
    unsigned long  biSize;
    long           biWidth;
    long           biHeight;
    unsigned short biPlanes;
    unsigned short biBitCount;
    unsigned long  biCompression;
    unsigned long  biSizeImage;
    long           biXPixPerMeter;
    long           biYPixPerMeter;
    unsigned long  biClrUsed;
    unsigned long  biClrImporant;
} BITMAPINFOHEADER;

#pragma pack(pop)

int exportToBmp(
	const char* fileName, 
	unsigned char* pixel, 
	unsigned int width,
	unsigned int height )
{
	//
	int success = 0;
	unsigned int x = 0;
	unsigned int y = 0;
	FILE* file = 0;
	BITMAPFILEHEADER* header = 0;
	BITMAPINFOHEADER* infoHeader = 0;
	unsigned int scanLineLengthInBytes = 0;
	int totalFileLength = 0;
	unsigned char* bmpMemoryStart = 0;
	unsigned char* bmpMemoryCursor = 0;
	//
	// スキャンライン長計算
	scanLineLengthInBytes = width*3;
	// ファイル全体の長さを計算
	totalFileLength = int(sizeof(BITMAPFILEHEADER)+sizeof(BITMAPINFOHEADER)+scanLineLengthInBytes*height);
	// メモリを確保
	bmpMemoryStart = (unsigned char*)malloc(totalFileLength);
	if( !bmpMemoryStart )
	{ goto EXIT; }
	bmpMemoryCursor = bmpMemoryStart;
	// BITMAPFILEHEADERを作成
	header = (BITMAPFILEHEADER*)bmpMemoryCursor; //
	header->bfType = 'B' | ('M' << 8); // ファイルタイプ
	header->bfSize = totalFileLength;// ファイルサイズ (byte)
	header->bfReserved1 = 0; // 予約領域
	header->bfReserved2 = 0; // 予約領域
	header->bfOffBits = sizeof(BITMAPFILEHEADER)+sizeof(BITMAPINFOHEADER); // ファイル先頭から画像データまでのオフセット (byte)
	bmpMemoryCursor += sizeof(BITMAPFILEHEADER);
	// BITMAPINFOHEADERを作成
	infoHeader = (BITMAPINFOHEADER*)bmpMemoryCursor;
	infoHeader->biSize = sizeof(BITMAPINFOHEADER); // 情報ヘッダのサイズ
	infoHeader->biWidth = width; // 画像の幅 (ピクセル)
	infoHeader->biHeight = height; // 画像の高さ (ピクセル)
	infoHeader->biPlanes = 1; // プレーン数
	infoHeader->biBitCount = 24; // 1 画素あたりのデータサイズ (bit)
	infoHeader->biCompression = 0; // 圧縮形式(無圧縮)
	infoHeader->biSizeImage = width*height*3; // 画像データ部のサイズ (byte)
	infoHeader->biXPixPerMeter = 3780; // 横方向解像度 (1mあたりの画素数)
	infoHeader->biYPixPerMeter = 3780; // 縦方向解像度 (1mあたりの画素数)
	infoHeader->biClrUsed = 0; // 格納されているパレット数 (使用色数)	
	infoHeader->biClrImporant = 0; // 重要なパレットのインデックス
	bmpMemoryCursor += sizeof(BITMAPINFOHEADER);
	// 全てのデータを書き込む
	for(y=0;y<height;++y)
	{
		for(x=0;x<width;++x)
		{
			const unsigned int srcBase = (x+(height-y-1)*width)*3;
			const unsigned int dstBase = x*3;
			bmpMemoryCursor[dstBase+0] = pixel[srcBase+0];
			bmpMemoryCursor[dstBase+1] = pixel[srcBase+1];
			bmpMemoryCursor[dstBase+2] = pixel[srcBase+2];
		}
		bmpMemoryCursor += scanLineLengthInBytes;
	}
	// ファイルに書き込む
#ifdef _WIN32
	fopen_s(&file,fileName,"wb");
#else
    file = fopen(fileName, "wb");
#endif
	if(!file)
	{ goto EXIT; }
	if(fwrite(bmpMemoryStart,totalFileLength,1,file)!=1)
	{ goto EXIT; }
	if(fclose(file)==EOF)
	{ goto EXIT; }
	// 最後まで到達したので成功
	success = 1;
EXIT:
	if(bmpMemoryStart)
	{
		free(bmpMemoryStart);
	}
	if(file)
	{
		fclose(file);
	}
	return success;
}
