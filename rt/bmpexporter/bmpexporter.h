#pragma once
#ifndef _BMP_EXPORTER_H_
#define _BMP_EXPORTER_H_

/**
24bitのピクセルデータを24bit形式のBMPに出力する
*/
extern int exportToBmp(
	const char* fileName, 
	unsigned char* pixel, 
	unsigned int width,
	unsigned int height );

#endif
