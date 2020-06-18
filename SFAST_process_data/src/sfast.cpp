#include <iostream>
#include "sfast.h"
//#include "hls_math.h"
#include "utils/x_hls_utils.h"

static col_pix_t glPLSlices[SLICES_NUMBER][SLICE_WIDTH][SLICE_HEIGHT/COMBINED_PIXELS];
static col_pix_t glPLSlicesScale1[SLICES_NUMBER][SLICE_WIDTH/2][SLICE_HEIGHT/COMBINED_PIXELS/2];
static col_pix_t glPLSlicesScale2[SLICES_NUMBER][SLICE_WIDTH/4][SLICE_HEIGHT/COMBINED_PIXELS/4];
static sliceIdx_t glPLActiveSliceIdx = 0, glPLTminus1SliceIdx, glPLTminus2SliceIdx;

static col_pix_t glPLSFASTSliceScale2[SLICES_NUMBER][SLICE_WIDTH/4][SLICE_HEIGHT/COMBINED_PIXELS/4];

//const int innerCircleOffset[INNER_SIZE][2] = {{0, 1}, {1, 1}, {1, 0}, {1, -1},
//		{0, -1}, {-1, -1}, {-1, 0}, {-1, 1}};
//const int outerCircleOffset[OUTER_SIZE][2] = {{0, 2}, {1, 2}, {2, 1}, {2, 0},
//		{2, -1}, {1, -2}, {0, -2}, {-1, -2},
//		{-2, -1}, {-2, 0}, {-2, 1}, {-1, 2}};

//const ap_int<64> innerTest =  ap_int<64>("0111101f0ffff0f1", 16);
//const ap_int<96> outerTest = ap_int<96>("021221202f1e0efeefe0e1f2", 16);

const ap_int<96> innerTest =  ap_int<96>("1f0ffff0f1011110", 16);
const ap_int<96> outerTest = ap_int<96>("2f1e0efeefe0e1f202122120", 16);

static ap_uint<2> glStage = 0;
static ap_uint<2> glStageBak = 0;
static hls::stream< ap_uint<2> >  glStageInStream("stageInStream");

static uint32_t glInitCounter = 0, glFeedbackCounter = 0;

// It takes 8192 cycles to reset the whole slice, so we use 13 bits.
static ap_uint<13> resetCnt;

pix_t readPixFromCol(col_pix_t colData, ap_uint<8> idx)
{
#pragma HLS INLINE
	pix_t retData;
	// Use bit selection plus for-loop to read multi-bits from a wider bit width value
	// rather than use range selection directly. The reason is that the latter will use
	// a lot of shift-register which will increase a lot of LUTs consumed.
	readWiderBitsLoop: for(int8_t yIndex = 0; yIndex < BITS_PER_PIXEL; yIndex++)
	{
#pragma HLS UNROLL
		const int bitOffset = BITS_PER_PIXEL >> 1;
		ap_uint<8 + bitOffset> colIdx;
		// Concatenate and bit shift rather than multiple and accumulation (MAC) can save area.
		colIdx.range(8 + bitOffset - 1, bitOffset) = ap_uint<10>(idx * BITS_PER_PIXEL).range(8 + bitOffset - 1, bitOffset);
		colIdx.range(bitOffset - 1, 0) = ap_uint<2>(yIndex);

		retData[yIndex] = colData[colIdx];
//		retData[yIndex] = colData[BITS_PER_PIXEL*idx + yIndex];
	}
	return retData;
}

pix_t readPixFromTwoCols(two_cols_pix_t colData, ap_uint<8> idx)
{
#pragma HLS INLINE
	pix_t retData;
	// Use bit selection plus for-loop to read multi-bits from a wider bit width value
	// rather than use range selection directly. The reason is that the latter will use
	// a lot of shift-register which will increase a lot of LUTs consumed.
//	ap_uint<256> colIdxHi, colIdxLo;
//	colIdxHi = (ap_uint<8>(idx * BITS_PER_PIXEL)(8,2), ap_uint<2>(0));
//	colIdxLo = (ap_uint<8>(idx * BITS_PER_PIXEL)(8,2), ap_uint<2>(BITS_PER_PIXEL - 1));
//	retData = colData(colIdxHi, colIdxLo);
	readTwoColsWiderBitsLoop: for(int8_t yIndex = 0; yIndex < BITS_PER_PIXEL; yIndex++)
	{
#pragma HLS UNROLL
		const int bitOffset = BITS_PER_PIXEL >> 1;
		ap_uint<8 + bitOffset> colIdx;
		// Concatenate and bit shift rather than multiple and accumulation (MAC) can save area.
		colIdx.range(8 + bitOffset - 1, bitOffset) = ap_uint<10>(idx * BITS_PER_PIXEL).range(8 + bitOffset - 1, bitOffset);
		colIdx.range(bitOffset - 1, 0) = ap_uint<2>(yIndex);

		retData[yIndex] = colData[colIdx];
//		retData[yIndex] = colData[BITS_PER_PIXEL*idx + yIndex];
	}
	return retData;
}

void writePixToCol(col_pix_t *colData, ap_uint<8> idx, pix_t pixData)
{
#pragma HLS INLINE
	writeWiderBitsLoop: for(int8_t yIndex = 0; yIndex < BITS_PER_PIXEL; yIndex++)
	{
#pragma HLS UNROLL
		const int bitOffset = BITS_PER_PIXEL >> 1;
		ap_uint<8 + bitOffset> colIdx;
		// Concatenate and bit shift rather than multiple and accumulation (MAC) can save area.
		colIdx.range(8 + bitOffset - 1, bitOffset) = ap_uint<10>(idx * BITS_PER_PIXEL).range(8 + bitOffset - 1, bitOffset);
		colIdx.range(bitOffset - 1, 0) = ap_uint<2>(yIndex);

		(*colData)[colIdx] = pixData[yIndex];
	}
}

void resetPix(apUint10_t x, apUint10_t y, sliceIdx_t sliceIdx)
{
#pragma HLS INLINE
	glPLSlices[sliceIdx][x][y/COMBINED_PIXELS] = 0;
	glPLSlicesScale1[sliceIdx][x/2][y/COMBINED_PIXELS/2] = 0;
	glPLSlicesScale2[sliceIdx][x/4][y/COMBINED_PIXELS/4] = 0;
	glPLSFASTSliceScale2[sliceIdx][x/4][y/COMBINED_PIXELS/4] = 0;
}


void writePix(apUint10_t x, apUint10_t y, sliceIdx_t sliceIdx)
{
#pragma HLS DEPENDENCE variable=glPLSlicesScale2 inter false
#pragma HLS DEPENDENCE variable=glPLSlicesScale1 inter false
#pragma HLS ARRAY_PARTITION variable=glPLSlicesScale2 cyclic factor=1 dim=3
#pragma HLS ARRAY_PARTITION variable=glPLSlicesScale2 complete dim=1
#pragma HLS RESOURCE variable=glPLSlicesScale1 core=RAM_T2P_BRAM
#pragma HLS ARRAY_PARTITION variable=glPLSlicesScale1 cyclic factor=1 dim=3
#pragma HLS ARRAY_PARTITION variable=glPLSlicesScale1 complete dim=1
#pragma HLS PIPELINE
#pragma HLS RESOURCE variable=glPLSlices core=RAM_T2P_BRAM
#pragma HLS INLINE
#pragma HLS ARRAY_PARTITION variable=glPLSlices complete dim=1
#pragma HLS ARRAY_PARTITION variable=glPLSlices cyclic factor=1 dim=3
#pragma HLS DEPENDENCE variable=glPLSlices inter false
	col_pix_t tmpData;
	pix_t tmpTmpData;

	ap_uint<8> yNewIdx = y%COMBINED_PIXELS;

	tmpData = glPLSlices[sliceIdx][x][y/COMBINED_PIXELS];

	tmpTmpData = readPixFromCol(tmpData, yNewIdx);

	ap_uint<1> cmpFlg = ap_uint<1>(tmpTmpData < (ap_uint< BITS_PER_PIXEL >(0xff)));
	tmpTmpData +=  cmpFlg;

	writePixToCol(&tmpData, yNewIdx, tmpTmpData);

	glPLSlices[sliceIdx][x][y/COMBINED_PIXELS] = tmpData;

    // write scale 1
	apUint10_t xScale1 = x/2;
	apUint10_t yScale1 = y/2;
    ap_uint<8> yNewIdxScale1 = yScale1%COMBINED_PIXELS;

	col_pix_t tmpDataScale1;
	pix_t tmpTmpDataScale1;

	tmpDataScale1 = glPLSlicesScale1[sliceIdx][xScale1][yScale1/COMBINED_PIXELS];
	tmpTmpDataScale1 = readPixFromCol(tmpDataScale1, yNewIdxScale1);
	ap_uint<1> cmpFlgScale1 = ap_uint<1>(tmpTmpDataScale1 < (ap_uint< BITS_PER_PIXEL >(0xff)));
	tmpTmpDataScale1 +=  cmpFlgScale1;
	writePixToCol(&tmpDataScale1, yNewIdxScale1, tmpTmpDataScale1);
    glPLSlicesScale1[sliceIdx][xScale1][yScale1/COMBINED_PIXELS] = tmpDataScale1;

    // write scale 2
    apUint10_t xScale2 = x/4;
    apUint10_t yScale2 = y/4;
    ap_uint<8> yNewIdxScale2 = yScale2%COMBINED_PIXELS;

	col_pix_t tmpDataScale2;
	pix_t tmpTmpDataScale2;

	tmpDataScale2 = glPLSlicesScale2[sliceIdx][xScale2][yScale2/COMBINED_PIXELS];
	tmpTmpDataScale2 = readPixFromCol(tmpDataScale2, yNewIdxScale2);
	ap_uint<1> cmpFlgScale2 = ap_uint<1>(tmpTmpDataScale2 < (ap_uint< BITS_PER_PIXEL >(0xff)));
	tmpTmpDataScale2 +=  cmpFlgScale2;
	writePixToCol(&tmpDataScale2, yNewIdxScale2, tmpTmpDataScale2);
    glPLSlicesScale2[sliceIdx][xScale2][yScale2/COMBINED_PIXELS] = tmpDataScale2;
}

void resetSFASTPix(apUint10_t x, apUint10_t y, sliceIdx_t sliceIdx)
{
#pragma HLS INLINE
	glPLSFASTSliceScale2[sliceIdx][y][x/COMBINED_PIXELS] = 0;
}

void writeSFASTPix(apUint10_t x, apUint10_t y, sliceIdx_t sliceIdx)
{
#pragma HLS DEPENDENCE variable=glPLSFASTSliceScale2 inter false
#pragma HLS RESOURCE variable=glPLSFASTSliceScale2 core=RAM_T2P_BRAM
#pragma HLS ARRAY_PARTITION variable=glPLSFASTSliceScale2 cyclic factor=1 dim=3
#pragma HLS ARRAY_PARTITION variable=glPLSFASTSliceScale2 complete dim=1
    // write sfast slice
    apUint10_t xSFASTScale2 = x/4;
    apUint10_t ySFASTScale2 = y/4;
    ap_uint<8> xNewIdxSFASTScale2 = xSFASTScale2%COMBINED_PIXELS;

	col_pix_t tmpDataSFASTScale2;
	pix_t tmpTmpDataSFASTScale2;

	tmpDataSFASTScale2 = glPLSFASTSliceScale2[sliceIdx][ySFASTScale2][xSFASTScale2/COMBINED_PIXELS];
	tmpTmpDataSFASTScale2 = readPixFromCol(tmpDataSFASTScale2, xNewIdxSFASTScale2);
	ap_uint<1> cmpFlgSFASTScale2 = ap_uint<1>(tmpTmpDataSFASTScale2 < (ap_uint< BITS_PER_PIXEL >(0xff)));
	tmpTmpDataSFASTScale2 +=  cmpFlgSFASTScale2;
	writePixToCol(&tmpDataSFASTScale2, xNewIdxSFASTScale2, tmpTmpDataSFASTScale2);
	glPLSFASTSliceScale2[sliceIdx][ySFASTScale2][xSFASTScale2/COMBINED_PIXELS] = tmpDataSFASTScale2;
}


void readBlockCols(apUint10_t x, apUint10_t y, sliceIdx_t sliceIdxRef, sliceIdx_t sliceIdxTag,
		pix_t refCol[BLOCK_SIZE + 2 * SEARCH_DISTANCE],
		pix_t tagCol[BLOCK_SIZE + 2 * SEARCH_DISTANCE])
{
#pragma HLS INLINE
#pragma HLS PIPELINE
#pragma HLS ARRAY_RESHAPE variable=refCol complete dim=1
#pragma HLS ARRAY_RESHAPE variable=tagCol complete dim=1

	two_cols_pix_t refColData;
    two_cols_pix_t tagColData;
    ap_uint<8> neighboryOffset;
    if ( y%COMBINED_PIXELS < BLOCK_SIZE/2 + SEARCH_DISTANCE )
    {
        neighboryOffset = y/COMBINED_PIXELS - 1;
    }
    else
    {
        neighboryOffset = y/COMBINED_PIXELS + 1;
    }

    // concatenate two columns together
    refColData = (glPLSlices[sliceIdxRef][x][y/COMBINED_PIXELS], glPLSlices[sliceIdxRef][x][neighboryOffset]);
    //	cout << "refColData: " << refColData << endl;
    // concatenate two columns together
    // Use explicit cast here, otherwise it will generate a lot of select operations which consumes more LUTs than MUXs.
    tagColData = (glPLSlices[(sliceIdx_t)(sliceIdxTag + 0)][x][y/COMBINED_PIXELS], glPLSlices[(sliceIdx_t)(sliceIdxTag + 0)][x][neighboryOffset]);

	// find the bottom pixel of the column that centered on y.
	ap_uint<6> yColOffsetIdx = y%COMBINED_PIXELS - BLOCK_SIZE/2 - SEARCH_DISTANCE + COMBINED_PIXELS;

	readRefLoop: for(ap_uint<8> i = 0; i < BLOCK_SIZE + 2 * SEARCH_DISTANCE; i++)
	{
		refCol[i] = readPixFromTwoCols(refColData,  yColOffsetIdx);
		tagCol[i] = readPixFromTwoCols(tagColData,  yColOffsetIdx);
		yColOffsetIdx++;
	}

}


void readBlockColsScale1(apUint10_t x, apUint10_t y, sliceIdx_t sliceIdxRef, sliceIdx_t sliceIdxTag,
		pix_t refColScale1[BLOCK_SIZE + 2 * SEARCH_DISTANCE],
		pix_t tagColScale1[BLOCK_SIZE + 2 * SEARCH_DISTANCE])
{
#pragma HLS ARRAY_RESHAPE variable=tagColScale1 complete dim=1
#pragma HLS ARRAY_RESHAPE variable=refColScale1 complete dim=1
#pragma HLS PIPELINE
#pragma HLS INLINE

	two_cols_pix_t refColData;
    two_cols_pix_t tagColData;
    ap_uint<8> neighboryOffset;
    if ( y%COMBINED_PIXELS < BLOCK_SIZE/2 + SEARCH_DISTANCE )
    {
        neighboryOffset = y/COMBINED_PIXELS - 1;
    }
    else
    {
        neighboryOffset = y/COMBINED_PIXELS + 1;
    }

    // concatenate two columns together
    refColData = (glPLSlicesScale1[sliceIdxRef][x][y/COMBINED_PIXELS], glPLSlicesScale1[sliceIdxRef][x][neighboryOffset]);
    //	cout << "refColData: " << refColData << endl;
    // concatenate two columns together
    // Use explicit cast here, otherwise it will generate a lot of select operations which consumes more LUTs than MUXs.
    tagColData = (glPLSlicesScale1[(sliceIdx_t)(sliceIdxTag + 0)][x][y/COMBINED_PIXELS], glPLSlicesScale1[(sliceIdx_t)(sliceIdxTag + 0)][x][neighboryOffset]);

	// find the bottom pixel of the column that centered on y.
	ap_uint<6> yColOffsetIdx = y%COMBINED_PIXELS - BLOCK_SIZE/2 - SEARCH_DISTANCE + COMBINED_PIXELS;

	readRefLoop: for(ap_uint<8> i = 0; i < BLOCK_SIZE + 2 * SEARCH_DISTANCE; i++)
	{
		refColScale1[i] = readPixFromTwoCols(refColData,  yColOffsetIdx);
		tagColScale1[i] = readPixFromTwoCols(tagColData,  yColOffsetIdx);
		yColOffsetIdx++;
	}
}

void readBlockColsScale2(apUint10_t x, apUint10_t y, sliceIdx_t sliceIdxRef, sliceIdx_t sliceIdxTag,
		pix_t refColScale2[BLOCK_SIZE + 2 * SEARCH_DISTANCE],
		pix_t tagColScale2[BLOCK_SIZE + 2 * SEARCH_DISTANCE])
{
#pragma HLS ARRAY_RESHAPE variable=tagColScale2 complete dim=1
#pragma HLS ARRAY_RESHAPE variable=refColScale2 complete dim=1
#pragma HLS PIPELINE
#pragma HLS INLINE

	two_cols_pix_t refColData;
    two_cols_pix_t tagColData;
    ap_uint<8> neighboryOffset;
    if ( y%COMBINED_PIXELS < BLOCK_SIZE/2 + SEARCH_DISTANCE )
    {
        neighboryOffset = y/COMBINED_PIXELS - 1;
    }
    else
    {
        neighboryOffset = y/COMBINED_PIXELS + 1;
    }

    // concatenate two columns together
    refColData = (glPLSlicesScale2[sliceIdxRef][x][y/COMBINED_PIXELS], glPLSlicesScale2[sliceIdxRef][x][neighboryOffset]);
    //	cout << "refColData: " << refColData << endl;
    // concatenate two columns together
    // Use explicit cast here, otherwise it will generate a lot of select operations which consumes more LUTs than MUXs.
    tagColData = (glPLSlicesScale2[(sliceIdx_t)(sliceIdxTag + 0)][x][y/COMBINED_PIXELS], glPLSlicesScale2[(sliceIdx_t)(sliceIdxTag + 0)][x][neighboryOffset]);

	// find the bottom pixel of the column that centered on y.
	ap_uint<6> yColOffsetIdx = y%COMBINED_PIXELS - BLOCK_SIZE/2 - SEARCH_DISTANCE + COMBINED_PIXELS;

	readRefLoop: for(ap_uint<8> i = 0; i < BLOCK_SIZE + 2 * SEARCH_DISTANCE; i++)
	{
		refColScale2[i] = readPixFromTwoCols(refColData,  yColOffsetIdx);
		tagColScale2[i] = readPixFromTwoCols(tagColData,  yColOffsetIdx);
		yColOffsetIdx++;
	}
}

void readBlockSFASTScale2(apUint10_t x, apUint10_t y, sliceIdx_t sliceIdxRef, sliceIdx_t sliceIdxTag,
		pix_t refColScale2[BLOCK_SIZE + 2 * SEARCH_DISTANCE],
		pix_t tagColScale2[BLOCK_SIZE + 2 * SEARCH_DISTANCE])
{
#pragma HLS ARRAY_RESHAPE variable=tagColScale2 complete dim=1
#pragma HLS ARRAY_RESHAPE variable=refColScale2 complete dim=1
#pragma HLS PIPELINE
#pragma HLS INLINE

	two_cols_pix_t refColData;
    two_cols_pix_t tagColData;
    ap_uint<8> neighboryOffset;
    if ( y%COMBINED_PIXELS < BLOCK_SIZE/2 + SEARCH_DISTANCE )
    {
        neighboryOffset = y/COMBINED_PIXELS - 1;
    }
    else
    {
        neighboryOffset = y/COMBINED_PIXELS + 1;
    }

    // concatenate two columns together
    refColData = (glPLSFASTSliceScale2[sliceIdxRef][x][y/COMBINED_PIXELS], glPLSFASTSliceScale2[sliceIdxRef][x][neighboryOffset]);
    //	cout << "refColData: " << refColData << endl;
    // concatenate two columns together
    // Use explicit cast here, otherwise it will generate a lot of select operations which consumes more LUTs than MUXs.
    tagColData = (glPLSFASTSliceScale2[(sliceIdx_t)(sliceIdxTag + 0)][x][y/COMBINED_PIXELS], glPLSFASTSliceScale2[(sliceIdx_t)(sliceIdxTag + 0)][x][neighboryOffset]);

	// find the bottom pixel of the column that centered on y.
	ap_uint<6> yColOffsetIdx = y%COMBINED_PIXELS - BLOCK_SIZE/2 - SEARCH_DISTANCE + COMBINED_PIXELS;

	readRefLoop: for(ap_uint<8> i = 0; i < BLOCK_SIZE + 2 * SEARCH_DISTANCE; i++)
	{
		refColScale2[i] = readPixFromTwoCols(refColData,  yColOffsetIdx);
		tagColScale2[i] = readPixFromTwoCols(tagColData,  yColOffsetIdx);
		yColOffsetIdx++;
	}
}

sliceIdx_t oldIdx = glPLActiveSliceIdx;
void rwSlices(hls::stream<apUint10_t> &xStream, hls::stream<apUint10_t> &yStream, hls::stream<sliceIdx_t> &idxStream,
			  hls::stream<apIntBlockCol_t> &refStreamOut, hls::stream<apIntBlockCol_t> &tagStreamOut,
			  hls::stream<apIntBlockCol_t> &refStreamOutScale1, hls::stream<apIntBlockCol_t> &tagStreamOutScale1,
			  hls::stream<apIntBlockCol_t> &refStreamOutScale2, hls::stream<apIntBlockCol_t> &tagStreamOutScale2)
{
#pragma HLS INLINE off
	apUint10_t xRd;
	apUint10_t yRd;
	sliceIdx_t idx;

	apIntBlockCol_t colData0[BLOCK_SIZE], colData1[BLOCK_SIZE + 2 * SEARCH_DISTANCE];

	pix_t out1_debug[BLOCK_SIZE + 2 * SEARCH_DISTANCE][BLOCK_SIZE + 2 * SEARCH_DISTANCE];
	pix_t out2_debug[BLOCK_SIZE + 2 * SEARCH_DISTANCE][BLOCK_SIZE + 2 * SEARCH_DISTANCE];

    pix_t out1Scale1_debug[BLOCK_SIZE + 2 * SEARCH_DISTANCE][BLOCK_SIZE+ 2 * SEARCH_DISTANCE];
    pix_t out2Scale1_debug[BLOCK_SIZE + 2 * SEARCH_DISTANCE][BLOCK_SIZE+ 2 * SEARCH_DISTANCE];

    pix_t out1Scale2_debug[BLOCK_SIZE + 2 * SEARCH_DISTANCE][BLOCK_SIZE+ 2 * SEARCH_DISTANCE];
    pix_t out2Scale2_debug[BLOCK_SIZE + 2 * SEARCH_DISTANCE][BLOCK_SIZE+ 2 * SEARCH_DISTANCE];

//	rwSlicesLoop:for(int32_t i = 0; i < eventIterSize; i++)
//	{
		rwSlicesInnerLoop:for(int8_t xOffSet = 0; xOffSet < BLOCK_SIZE * (2 * SEARCH_DISTANCE + 1); xOffSet++)
		{
#pragma HLS PIPELINE rewind
//			xRd = (xOffSet == 0)? (ap_uint<8>)(xStream.read()): xRd;
//			yRd = (xOffSet == 0)? (ap_uint<8>)(yStream.read()): yRd;
			if (xOffSet == 0)
			{
				xRd = xStream.read();
				yRd = yStream.read();
				idx = idxStream.read();

				/* This is only for C-simulation and debugging. */
				if (oldIdx != idx)
				{
					oldIdx = idx;
					// Check the accumulation slice is clear or not
					for(int32_t xAddr = 0; xAddr < SLICE_WIDTH; xAddr++)
					{
						for(int32_t yAddr = 0; yAddr < SLICE_HEIGHT; yAddr = yAddr + COMBINED_PIXELS)
						{
							if (glPLSlices[idx][xAddr][yAddr/COMBINED_PIXELS] != 0)
							{
								for(int r = 0; r < 1000; r++)
								{
									std::cout << "Ha! I caught you, the pixel which is not clear!" << std::endl;
									std::cout << "x is: " << xAddr << "\t y is: " << yAddr << "\t idx is: " << idx << std::endl;
								}
							}
						}
					}
				}

#ifdef DEBUG
				/* This is only for C-simulation and debugging. */
				debugRegion:
				{
					for(int i = 0; i < BLOCK_SIZE + 2 * SEARCH_DISTANCE; i++)
					{
						readBlockCols(xRd - BLOCK_SIZE/2 - SEARCH_DISTANCE + i, yRd , idx + 1, idx + 2, out1_debug[i], out2_debug[i]);
						readBlockColsScale1(xRd/2 - BLOCK_SIZE/2 - SEARCH_DISTANCE + i, yRd/2 , idx + 1, idx + 2, out1Scale1_debug[i], out2Scale1_debug[i]);
						readBlockColsScale2(xRd/4 - BLOCK_SIZE/2 - SEARCH_DISTANCE + i, yRd/4 , idx + 1, idx + 2, out1Scale2_debug[i], out2Scale2_debug[i]);
					}
				}
#endif

				writePix(xRd, yRd, idx);

				resetPix(resetCnt/(PIXS_PER_COL), (resetCnt % (PIXS_PER_COL)) * COMBINED_PIXELS, (sliceIdx_t)(idx + 3));
				resetCnt++;
			}
			else if(xOffSet < BLOCK_SIZE + 2 * SEARCH_DISTANCE + 1)
			{
				pix_t out1[BLOCK_SIZE + 2 * SEARCH_DISTANCE];
				pix_t out2[BLOCK_SIZE + 2 * SEARCH_DISTANCE];

	            pix_t out1Scale1[BLOCK_SIZE+ 2 * SEARCH_DISTANCE];
	            pix_t out2Scale1[BLOCK_SIZE+ 2 * SEARCH_DISTANCE];

	            pix_t out1Scale2[BLOCK_SIZE+ 2 * SEARCH_DISTANCE];
	            pix_t out2Scale2[BLOCK_SIZE+ 2 * SEARCH_DISTANCE];

				readBlockCols(xRd - BLOCK_SIZE/2 - SEARCH_DISTANCE + xOffSet - 1, yRd , idx + 1, idx + 2, out1, out2);
				readBlockColsScale1(xRd/2 - BLOCK_SIZE/2 - SEARCH_DISTANCE + xOffSet - 1, yRd/2 , idx + 1, idx + 2, out1Scale1, out2Scale1);
				readBlockColsScale2(xRd/4 - BLOCK_SIZE/2 - SEARCH_DISTANCE + xOffSet - 1, yRd/4 , idx + 1, idx + 2, out1Scale2, out2Scale2);

				apIntBlockCol_t refBlockCol;
				apIntBlockCol_t tagBlockCol;

				apIntBlockCol_t refBlockColScale1;
				apIntBlockCol_t tagBlockColScale1;

				apIntBlockCol_t refBlockColScale2;
				apIntBlockCol_t tagBlockColScale2;

				for (int8_t l = 0; l < BLOCK_SIZE + 2 * SEARCH_DISTANCE; l++)
				{
					refBlockCol.range(BITS_PER_PIXEL * l + BITS_PER_PIXEL - 1, BITS_PER_PIXEL * l) = out1[l];
					tagBlockCol.range(BITS_PER_PIXEL * l + BITS_PER_PIXEL - 1, BITS_PER_PIXEL * l) = out2[l];

					refBlockColScale1.range(BITS_PER_PIXEL * l + BITS_PER_PIXEL - 1, BITS_PER_PIXEL * l) = out1Scale1[l];
					tagBlockColScale1.range(BITS_PER_PIXEL * l + BITS_PER_PIXEL - 1, BITS_PER_PIXEL * l) = out2Scale1[l];

					refBlockColScale2.range(BITS_PER_PIXEL * l + BITS_PER_PIXEL - 1, BITS_PER_PIXEL * l) = out1Scale2[l];
					tagBlockColScale2.range(BITS_PER_PIXEL * l + BITS_PER_PIXEL - 1, BITS_PER_PIXEL * l) = out2Scale2[l];
				}

				if (xOffSet > SEARCH_DISTANCE && xOffSet <= SEARCH_DISTANCE + BLOCK_SIZE)
				{
					refStreamOut << refBlockCol;
					refStreamOutScale1 <<  refBlockColScale1;
					refStreamOutScale2 <<  refBlockColScale2;
				}
				tagStreamOut << tagBlockCol;
				tagStreamOutScale1 << tagBlockColScale1;
				tagStreamOutScale2 << tagBlockColScale2;
			}
			else
			{
				// Reset two pixels at the same time because it has two write ports.
				resetPix(resetCnt/(PIXS_PER_COL), (resetCnt % (PIXS_PER_COL)) * COMBINED_PIXELS, (sliceIdx_t)(idx + 3));
				resetCnt++;
				resetPix(resetCnt/(PIXS_PER_COL), (resetCnt % (PIXS_PER_COL)) * COMBINED_PIXELS, (sliceIdx_t)(idx + 3));
				resetCnt++;
			}
		}
//	}

#ifdef DEBUG
	printRegion: if(xRd == 211 && yRd == 242)
	{
		std::cout << "Reference block of scale 0 is: " << std::endl;
		for(uint8_t blockX = SEARCH_DISTANCE; blockX < BLOCK_SIZE + SEARCH_DISTANCE; blockX++)
		{
			for(uint8_t blockY = SEARCH_DISTANCE; blockY < BLOCK_SIZE + SEARCH_DISTANCE; blockY++)
			{
				std::cout << out1_debug[blockX][blockY] << "\t";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;

		std::cout << "target block of scale 0 is: " << std::endl;
		for(uint8_t blockX = 0; blockX < BLOCK_SIZE + 2 * SEARCH_DISTANCE; blockX++)
		{
			for(uint8_t blockY = 0; blockY < BLOCK_SIZE + 2 * SEARCH_DISTANCE; blockY++)
			{
				std::cout << out2_debug[blockX][blockY] << "\t";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;

		std::cout << "Reference block of scale 1 is: " << std::endl;
		for(uint8_t blockX = SEARCH_DISTANCE; blockX < BLOCK_SIZE + SEARCH_DISTANCE; blockX++)
		{
			for(uint8_t blockY = SEARCH_DISTANCE; blockY < BLOCK_SIZE + SEARCH_DISTANCE; blockY++)
			{
				std::cout << out1Scale1_debug[blockX][blockY] << "\t";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;

		std::cout << "target block of scale 1 is: " << std::endl;
		for(uint8_t blockX = 0; blockX < BLOCK_SIZE + 2 * SEARCH_DISTANCE; blockX++)
		{
			for(uint8_t blockY = 0; blockY < BLOCK_SIZE + 2 * SEARCH_DISTANCE; blockY++)
			{
				std::cout << out2Scale1_debug[blockX][blockY] << "\t";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;

		std::cout << "Reference block of scale 2 is: " << std::endl;
		for(uint8_t blockX = SEARCH_DISTANCE; blockX < BLOCK_SIZE + SEARCH_DISTANCE; blockX++)
		{
			for(uint8_t blockY = SEARCH_DISTANCE; blockY < BLOCK_SIZE + SEARCH_DISTANCE; blockY++)
			{
				std::cout << out1Scale2_debug[blockX][blockY] << "\t";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;

		std::cout << "target block of scale 2 is: " << std::endl;
		for(uint8_t blockX = 0; blockX < BLOCK_SIZE + 2 * SEARCH_DISTANCE; blockX++)
		{
			for(uint8_t blockY = 0; blockY < BLOCK_SIZE + 2 * SEARCH_DISTANCE; blockY++)
			{
				std::cout << out2Scale2_debug[blockX][blockY] << "\t";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}


	resetLoop: for (int16_t resetCnt = 0; resetCnt < SLICE_HEIGHT * SLICE_WIDTH / COMBINED_PIXELS; resetCnt = resetCnt + 2)
	{
		resetPix(resetCnt/PIXS_PER_COL, (resetCnt % PIXS_PER_COL) * COMBINED_PIXELS, (sliceIdx_t)(idx + 3));
		resetPix(resetCnt/PIXS_PER_COL, (resetCnt % PIXS_PER_COL + 1) * COMBINED_PIXELS, (sliceIdx_t)(idx + 3));
	}
#endif


//	resetLoop: for (int16_t resetCnt = 0; resetCnt < SLICE_HEIGHT * SLICE_WIDTH / COMBINED_PIXELS; resetCnt = resetCnt + 2)
//	{
//		resetPix(resetCnt/PIXS_PER_COL, (resetCnt % PIXS_PER_COL) * COMBINED_PIXELS, (sliceIdx_t)(idx + 3));
//		resetPix(resetCnt/PIXS_PER_COL, (resetCnt % PIXS_PER_COL + 1) * COMBINED_PIXELS, (sliceIdx_t)(idx + 3));
//	}

}

static ap_uint<32> glConfig;
static status_t glStatus;
static uint64_t inEventsNum;
void preProcessStream(hls::stream< ap_uint<16> > &xStreamIn, hls::stream< ap_uint<16> > &yStreamIn,
		hls::stream< ap_uint<1> > &polStreamIn, hls::stream<sliceIdx_t> &idxStreamIn,
		hls::stream< ap_uint<64> > &tsStreamIn,
		hls::stream<X_TYPE> &xStreamOut, hls::stream<Y_TYPE> &yStreamOut, hls::stream< ap_uint<1> > &polStreamOut,
		hls::stream<sliceIdx_t> &idxStreamOut,
		hls::stream< ap_uint<TS_TYPE_BIT_WIDTH> > &tsStreamOut, hls::stream< ap_uint<96> > &packetEventDataStream)
{
#pragma HLS PIPELINE
	ap_uint<16> x;
	ap_uint<16> y;
	ap_uint<64> ts;
	ap_uint<1> pol;

	xStreamIn >> x;
	yStreamIn >> y;
	tsStreamIn >> ts;
	polStreamIn >> pol;

	// A register to store current slice index and init it as 0
	// Update it only when the input idxStreamIn is not empty;
	static sliceIdx_t sliceIdxReg = 0;
	if(!idxStreamIn.empty())
	{
		sliceIdxReg = idxStreamIn.read();
	}
	else
	{
//		sliceIdxReg = sliceIdxReg;
	}

	ap_uint<96> tmpOutput;
	tmpOutput[32] = ap_uint<1>(pol);
	tmpOutput.range(31, 16) = y;
	tmpOutput.range(15, 0) = x;
	tmpOutput.range(95, 33) = ts.range(62, 0);
	packetEventDataStream << tmpOutput;

	// TODO: Removed the hardcoded code invalid event
	const int max_scale = 2;
	// only check if not too close to border
	const int cs = 4;
	// Make this event an invalid event
	// 4 is the corner block range.
	// Remember x,y is inverted in the raw stream.
	if ((x >> max_scale) < cs || (x >> max_scale) >= (DVS_REAL_WIDTH >> max_scale) - cs ||
			(y >> max_scale) < cs || (y >> max_scale) >= (DVS_REAL_HEIGHT >> max_scale) - cs)
	{
		x = 0;
		y = 0;
		ts = 0;
	}
	else
	{
		 // The reason why we sub cs-4 is to reduce the memory we need.
		 // In this case, only memory for valid events is required.
//		x -= cs-4;
//		y -= cs-4;
	}


	xStreamOut << (X_TYPE)x;
	yStreamOut << (Y_TYPE)y;
	polStreamOut << pol;
	tsStreamOut << (ap_uint<TS_TYPE_BIT_WIDTH>)ts;
	idxStreamOut << sliceIdxReg;

	//	xStreamOut << (X_TYPE)x;
//	yStreamOut << (Y_TYPE)y;
//	tsStreamOut << (ap_uint<TS_TYPE_BIT_WIDTH>)ts;

	inEventsNum++;
	glStatus.inEventsNum = inEventsNum;
}

void initStageInterleaveStream(ap_uint<2> *stageOutStream)
{
#pragma HLS INLINE off
	ap_uint<2> stageIn = 0;
	ap_uint<1> isStageCorner = 0;
	if(glInitCounter%(GROUP_EVENTS_NUM*2) < GROUP_EVENTS_NUM)
	{
		stageIn = 0;
	}
    else
    {
    	stageIn = glStageInStream.read();
//    	stageIn = 1;
    }
//	glStage = stageIn;

	// This counter is used to synchronize this module and the feedback module.
	// Make them have the common adder source.
	glFeedbackCounter = glInitCounter;
    glInitCounter++;

	*stageOutStream = stageIn;
}


template<int READ_NPC>   //  Due to the memory has 2 ports at most for arbitrary reading, here this number could be only 1 or 2.
void rwSAEPerfectLoopStreamV2(hls::stream<X_TYPE> &xStream, hls::stream<Y_TYPE> &yStream,
		hls::stream< ap_uint<1> > &polStream, hls::stream<sliceIdx_t> &idxStream,
		hls::stream< ap_uint<TS_TYPE_BIT_WIDTH> > &tsStream,
//		ap_uint<2> stageStream,
		hls::stream< ap_uint<BITS_PER_PIXEL * OUTER_SIZE> > &outputDataStream,
		ap_uint<5> *sizeStreamOut1,
		ap_uint<5> *sizeStreamOut2)
{
#pragma HLS INLINE off

	X_TYPE x;
	Y_TYPE y;
	ap_uint<1> pol;
	ap_uint<TS_TYPE_BIT_WIDTH> ts;
	ap_uint<2> stage = 1;
	sliceIdx_t idx;

	ap_uint<BITS_PER_PIXEL * OUTER_SIZE> tmpData = 0;
	ap_uint<5> size = OUTER_SIZE;

//	x = xStream.read();
//	y = yStream.read();
//	pol = polStream.read();
//	ts = tsStream.read();
//	idx = idxStream.read();

//	stage = stageStream.read();

	for(ap_int<8> i = -1 * READ_NPC; i < OUTER_SIZE; i = i + READ_NPC)
	{
#pragma HLS PIPELINE rewind

		if(i == -1 * READ_NPC)
		{
			x = xStream.read();
			y = yStream.read();
			pol = polStream.read();
			ts = tsStream.read();
			idx = idxStream.read();
			writeSFASTPix(x, y, idx);

			resetSFASTPix(resetCnt/(PIXS_PER_COL/4), (resetCnt % (PIXS_PER_COL/4)) * COMBINED_PIXELS, (sliceIdx_t)(idx + 3));
			resetCnt++;
		}
		else
		{
            ap_uint<8 * READ_NPC> xInnerTest, xOuterTest, indexUnified;
            X_TYPE xInner[READ_NPC], xInnerNewIdx[READ_NPC];
            Y_TYPE yInner[READ_NPC];
        	rwSAEReadOffsetBitsLoop:
            for (ap_uint<8> j = 0; j < 8 * READ_NPC; j++)
            {
#pragma HLS UNROLL
            	ap_uint<8> tmpIndex;   // In  order to save the resource, we use bit operation to get the index of the inner/outer offset.
            	tmpIndex.range(7, 2 + READ_NPC) = ap_uint<8>(i * 8).range(7, 2 + READ_NPC);
            	tmpIndex.range(1 + READ_NPC, 0) = j.range(1 + READ_NPC, 0);
            	xInnerTest[j] = innerTest[tmpIndex];
            	xOuterTest[j] = outerTest[tmpIndex];
            }

            if(stage == 0)
            {
            	indexUnified = xInnerTest;
            }
            else if(stage == 1)
            {
            	indexUnified = xOuterTest;
            }
//		            indexUnified = (stage==0) ? xInnerTest : xOuterTest;

            readNPCLoop:
            for (ap_uint<8> k = 0; k < READ_NPC; k++)
            {
            	pix_t tmpTmpData, tmpTmpData0, tmpTmpData1;

                xInner[k] = x/4 + ap_int<4>(indexUnified(8 * k + 3, 8  * k));  // Change back from unsigned to signed.
                yInner[k] = y/4 + ap_int<4>(indexUnified(8 * k + 7, 8 * k + 4));          // Change back from unsigned to signed.
                xInnerNewIdx[k] = xInner[k]%COMBINED_PIXELS;

                tmpTmpData = readPixFromCol(glPLSFASTSliceScale2[sliceIdx_t(idx + 1)][yInner[k]][xInner[k]/COMBINED_PIXELS], xInnerNewIdx[k]);

    			// By doing this, we can avoid swapping the data order.
                tmpData = tmpData >> BITS_PER_PIXEL;
                tmpData.range(BITS_PER_PIXEL * OUTER_SIZE - 1, BITS_PER_PIXEL * OUTER_SIZE - BITS_PER_PIXEL) = ((stage == 0) && (i >= INNER_SIZE)) ?  pix_t(0) : tmpTmpData;
            }

//			resetSFASTPix(resetCnt/(PIXS_PER_COL), (resetCnt % (PIXS_PER_COL)) * COMBINED_PIXELS, (sliceIdx_t)(idx + 3));
//			resetCnt++;
//			resetSFASTPix(resetCnt/(PIXS_PER_COL), (resetCnt % (PIXS_PER_COL)) * COMBINED_PIXELS, (sliceIdx_t)(idx + 3));
//			resetCnt++;
		}
	}

    outputDataStream << tmpData;

	*sizeStreamOut1 = size;
	*sizeStreamOut2 = size;
}

void testRWSAEPerfectLoopStreamV2(hls::stream<X_TYPE> &xStream, hls::stream<Y_TYPE> &yStream,
		hls::stream< ap_uint<1> > &polStream, hls::stream<sliceIdx_t> &idxStream,
		hls::stream< ap_uint<TS_TYPE_BIT_WIDTH> > &tsStream,
		hls::stream< ap_uint<BITS_PER_PIXEL * OUTER_SIZE> > &outputDataStream,
		ap_uint<5> *sizeStreamOut1,
		ap_uint<5> *sizeStreamOut2)
{
	rwSAEPerfectLoopStreamV2<2>(xStream, yStream, polStream, idxStream, tsStream, outputDataStream, sizeStreamOut1, sizeStreamOut2);
}

// Function Description: convert the current data array to sorted idx array.
// Compared with V2, this version use a lot of shift logic tricks and it saves a lot of LUTs (most of them are MUXs).
template<int NPC>
void sortedIdxStreamV3(hls::stream< ap_uint<BITS_PER_PIXEL * OUTER_SIZE> > &tsStream, ap_uint<5> size, ap_uint<5*OUTER_SIZE> *newIdx)
{
assert(size <= OUTER_SIZE);
#pragma HLS INLINE off
	ap_uint<BITS_PER_PIXEL * OUTER_SIZE> tmpData = tsStream.read();
	ap_uint<BITS_PER_PIXEL> inData[OUTER_SIZE];

	ap_uint<5*OUTER_SIZE> tmpDataRet = 0;
	ap_uint<BITS_PER_PIXEL * OUTER_SIZE> dataToBeComp = tmpData;

	for(uint8_t i = 0; i < OUTER_SIZE/NPC; i = i + 1)
	{
// #pragma HLS LOOP_TRIPCOUNT min=0 max=20/NPC
#pragma HLS PIPELINE rewind
		if (0 == size)
		{
			break;
		}
		for(uint8_t j = 0; j < NPC; j++)
		{
			ap_uint<5> tmpIdx = 0;
			ap_uint<BITS_PER_PIXEL * OUTER_SIZE> tsWideData = tmpData;

			// We always shift the data to 5 MSBs and then shift right.
			// By doing this, we can avoid swapping the data order.
			tmpDataRet = tmpDataRet >> 5;
//			idxSorted(dataToBeComp.range(TS_TYPE_BIT_WIDTH - 1, 0), inData, size, &tmpIdx);

			sortedRegion:
			{
				for(uint8_t i = 0; i < OUTER_SIZE; i++ )
				{
					ap_uint<1> cond1 = (tsWideData.range(BITS_PER_PIXEL - 1, 0) < dataToBeComp.range(BITS_PER_PIXEL - 1, 0));  // Notice the difference between < and <= here.
					tmpIdx += cond1;
					tsWideData = tsWideData >> BITS_PER_PIXEL;
				}
			}
			dataToBeComp = dataToBeComp >> BITS_PER_PIXEL;
			tmpDataRet.range(5*OUTER_SIZE - 1, 5*OUTER_SIZE - 5) = tmpIdx;
			*newIdx = tmpDataRet;
		}
	}
}

template<int NPC>
void checkIdxGeneralV3(ap_uint<5*OUTER_SIZE> idxData, ap_uint<5> size, ap_uint<1> *isCorner)
{
#pragma HLS INLINE off
#pragma HLS ARRAY_PARTITION variable=idxData complete dim=0

	ap_uint<1> isCornerTemp = 0;
	ap_uint<5*INNER_SIZE> tmpIdxInnerData;
	ap_uint<5*OUTER_SIZE> tmpIdxOuterData;
	ap_uint<5> idxInnerData[INNER_SIZE];
	ap_uint<5> idxOuterData[OUTER_SIZE];

	for(uint8_t i = 0; i < OUTER_SIZE; i = i + NPC)
	{
#pragma HLS PIPELINE rewind
		InitRegion:
		{
			if (i == 0)
			{
				tmpIdxInnerData = idxData;
				tmpIdxOuterData = idxData;
			}
		}

		if(size == INNER_SIZE)
		{
			ap_uint<1> cond[4][6 + NPC - 1];
			for (uint8_t m = 0; m < 3 + NPC - 1; m++)
			{
				// The condition should be the idxData > (INNER_SIZE -3).
				// However, in order to make the idxSorted could be shared by inner circle and outer circle together.
				// We use a method that compare "size" values to all the input data which has OUTER_SIZE values in total.
				// On the other hand, if the valid input data number is less than OUTER_SIZE, the other input data will be filled with 0.
				// Thus, all the idxData for inner circle value will be added 4 (OUTER_SIZE - INNER_SIZE = 20 - 16 =4)
				// When we check the innner idx data, we need to remove it.
				cond[0][m] = (tmpIdxInnerData.range(5 * m + 4, 5 * m) >= INNER_SIZE - 3 + OUTER_SIZE - INNER_SIZE);
			}

			ap_uint<1> cond2[4 + NPC - 1];
			for (uint8_t m = 0; m < 4 + NPC - 1; m++)
			{
				cond[1][m] = (tmpIdxInnerData.range(5 * m + 4, 5 * m) >= INNER_SIZE - 4 + OUTER_SIZE - INNER_SIZE);
			}

			ap_uint<1> cond3[5 + NPC - 1];
			for (uint8_t m = 0; m < 5 + NPC - 1; m++)
			{
				cond[2][m] = (tmpIdxInnerData.range(5 * m + 4, 5 * m) >= INNER_SIZE - 5 + OUTER_SIZE - INNER_SIZE);
			}

			ap_uint<1> cond4[6 + NPC - 1];
			for (uint8_t m = 0; m < 6 + NPC - 1; m++)
			{
				cond[3][m] = (tmpIdxInnerData.range(5 * m + 4, 5 * m) >= INNER_SIZE - 6 + OUTER_SIZE - INNER_SIZE);
			}

			ap_uint<1> tempCond[4][NPC];

			for (uint8_t k = 0; k < NPC; k++)
			{
				for (uint8_t n = 0; n < 4; n++)
				{
					tempCond[n][k] = 1;
					for (uint8_t j = 0; j < 3 + n; j++)
					{
						tempCond[n][k] &= cond[n][j + k];
					}
					isCornerTemp |= tempCond[n][k];
	//				if (isCornerTemp == 1)
	//				{
	//					*isCorner = isCornerTemp ;
	//					std::cout << "HW: Position is :" << (int)(i + k) << " and streak size is: " << (int)(n + 3) << std::endl;
	//					return;
	//				}
				}
			}

			// Cycle shift the data list NPC times since we have processed NPC data.
			// In this loop, every loop shift the data one time.
			for(int shiftCnt = 0; shiftCnt < NPC; shiftCnt++)
			{
				tmpIdxInnerData = (tmpIdxInnerData >> 5) | (tmpIdxInnerData << (5*INNER_SIZE - 5));
			}

		}
		else if(size == OUTER_SIZE)
		{
			ap_uint<1> cond[5][8 + NPC - 1];
			for (uint8_t m = 0; m < 4 + NPC - 1; m++)
			{
				cond[0][m] = (tmpIdxOuterData.range(5 * m + 4, 5 * m) >= OUTER_SIZE - 4);
			}

			ap_uint<1> cond2[5 + NPC - 1];
			for (uint8_t m = 0; m < 5 + NPC - 1; m++)
			{
				cond[1][m] = (tmpIdxOuterData.range(5 * m + 4, 5 * m) >= OUTER_SIZE - 5);
			}

			ap_uint<1> cond3[6 + NPC - 1];
			for (uint8_t m = 0; m < 6 + NPC - 1; m++)
			{
				cond[2][m] = (tmpIdxOuterData.range(5 * m + 4, 5 * m) >= OUTER_SIZE - 6);
			}

			ap_uint<1> cond4[7 + NPC - 1];
			for (uint8_t m = 0; m < 7 + NPC - 1; m++)
			{
				cond[3][m] = (tmpIdxOuterData.range(5 * m + 4, 5 * m) >= OUTER_SIZE - 7);
			}

			ap_uint<1> cond5[8 + NPC - 1];
			for (uint8_t m = 0; m < 8 + NPC - 1; m++)
			{
				cond[4][m] = (tmpIdxOuterData.range(5 * m + 4, 5 * m) >= OUTER_SIZE - 8);
			}

			ap_uint<1> tempCond[5][NPC];

			for (uint8_t k = 0; k < NPC; k++)
			{
				for (uint8_t n = 0; n < 5; n++)
				{
					tempCond[n][k] = 1;
					for (uint8_t j = 0; j < 4 + n; j++)
					{
						tempCond[n][k] &= cond[n][j + k];
					}
					isCornerTemp |= tempCond[n][k];
	//				if (isCornerTemp == 1)
	//				{
	//					*isCorner = isCornerTemp ;
	//					std::cout << "HW: Position is :" << (int)(i + k) << " and streak size is: " << (int)(n + 4) << std::endl;
	//					return;
	//				}
				}
			}

			// Cycle shift the data list NPC times since we have processed NPC data.
			// In this loop, every loop shift the data one time.
			for(int shiftCnt = 0; shiftCnt < NPC; shiftCnt++)
			{
				tmpIdxOuterData = (tmpIdxOuterData >> 5) | (tmpIdxOuterData << (5*OUTER_SIZE - 5));
			}
		}
		else
		{
			break;
		}
	}
	*isCorner = isCornerTemp;
}

void feedbackInterleaveStream(ap_uint<1> isStageCorner, hls::stream< ap_uint<1> > &isFinalCornerStream)
{
#pragma HLS INLINE off
#pragma PIPELINE rewind
	ap_uint<2> outputStage;

	static ap_uint<2> currentStage[GROUP_EVENTS_NUM];
#pragma HLS ARRAY_PARTITION variable=currentStage complete dim=1

	ap_uint<2> tmpStage = currentStage[glFeedbackCounter%GROUP_EVENTS_NUM];

	if(tmpStage == 0)
	{
		if(isStageCorner)
		{
			tmpStage =  1;
		}
		else
		{
			tmpStage =  2;
		}
	}
	else
	{
		tmpStage =  0;
	}


	if(glFeedbackCounter%(GROUP_EVENTS_NUM * 2) < GROUP_EVENTS_NUM)
	{
//		stageStream << outputStage;
//		glStageInStream << tmpStage;
		isFinalCornerStream << isStageCorner;
	}
	else
	{
//		glStageInStream << glStage;
		isFinalCornerStream << isStageCorner;
	}
}

static uint64_t outEventsNum, cornerEventsNum;

void combineOutputStream(hls::stream< ap_uint<96> > &packetEventDataStream, hls::stream< ap_uint<1> > &stageCornerStream,
						hls::stream< ap_uint<16> > &xStreamOut, hls::stream< ap_uint<16> > &yStreamOut,
						hls::stream< ap_uint<1> > &polStreamOut,
						hls::stream< ap_uint<64> > &tsStreamOut, hls::stream< ap_uint<1> > &isFinalCornerStream)
{
#pragma HLS PIPELINE
	ap_uint<96> tmpOutput;
	packetEventDataStream >> tmpOutput;
	ap_uint<16> x;
	ap_uint<16> y;
	ap_uint<64> ts;
	ap_uint<1> pol;

	y = tmpOutput.range(31, 16);
	x = tmpOutput.range(15, 0);
	pol = tmpOutput[32];
	ts.range(62, 0) = tmpOutput.range(95, 33);

//	ap_uint<1> isCornerStage0 = isFinalCornerStream.read();
//	ap_uint<1> isCornerStage1 = isFinalCornerStream.read();
//	ap_uint<1> cornerRet = isCornerStage0 & isCornerStage1;

	ap_uint<1> cornerRet = stageCornerStream.read();

	ap_uint<8> pixelData;
	pixelData = (cornerRet == 1) ? 0xaa : 0;

	ap_uint<10> retData;
	retData = cornerRet;

//	if(glConfig[0])                   // This is filter mode.
//	{
//		if(cornerRet == 1 || ts.range(14,0) == 0) // time wrapping events should not be skipped
//		{
//			xStreamOut << x;
//			yStreamOut << y;
//			polStreamOut << pol;
//			tsStreamOut << ts;
//			custDataStreamOut << pixelData;
//
//		}
//	}
//	else                           // This is forward mode.
//	{
//		xStreamOut << x;
//		yStreamOut << y;
//		polStreamOut << pol;
//		tsStreamOut << ts;
//		custDataStreamOut << pixelData;
//	}

	xStreamOut << x;
	yStreamOut << y;
	polStreamOut << pol;
	tsStreamOut << ts;
	isFinalCornerStream << cornerRet;

	if(cornerRet == 1)
	{
		cornerEventsNum++;
	}
//	outEventsNum++;
//	glStatus.cornerEventsNum = cornerEventsNum;
//	glStatus.outEventsNum = outEventsNum;
}


void SFAST_process_data(hls::stream< ap_uint<16> > &xStreamIn, hls::stream< ap_uint<16> > &yStreamIn,
		hls::stream< ap_uint<64> > &tsStreamIn, hls::stream< ap_uint<1> > &polStreamIn, hls::stream<sliceIdx_t> &idxStreamIn,
		hls::stream< ap_uint<16> > &xStreamOut, hls::stream< ap_uint<16> > &yStreamOut, hls::stream< ap_uint<64> > &tsStreamOut, hls::stream< ap_uint<1> > &polStreamOut,
		hls::stream< ap_uint<1> > &isFinalCornerStream)
//		ap_uint<32> config, status_t *status)
{
#pragma HLS STREAM variable=glStageInStream depth=10 dim=1
#pragma HLS INTERFACE axis register both port=isFinalCornerStream
//#pragma HLS INTERFACE s_axilite port=status bundle=config
//#pragma HLS INTERFACE s_axilite port=config bundle=config
#pragma HLS INTERFACE axis register both port=tsStreamOut
#pragma HLS INTERFACE axis register both port=polStreamOut
#pragma HLS INTERFACE axis register both port=yStreamOut
#pragma HLS INTERFACE axis register both port=xStreamOut

#pragma HLS INTERFACE axis register both port=tsStreamIn
#pragma HLS INTERFACE axis register both port=yStreamIn
#pragma HLS INTERFACE axis register both port=xStreamIn
#pragma HLS INTERFACE axis register both port=polStreamIn
#pragma HLS DATAFLOW

	ap_uint<TS_TYPE_BIT_WIDTH> outer[OUTER_SIZE];

	hls::stream< ap_uint<BITS_PER_PIXEL * OUTER_SIZE> > inStream("dataStream");
#pragma HLS STREAM variable=inStream depth=10 dim=1
#pragma HLS RESOURCE variable=inStream core=FIFO_SRL

	ap_uint<5> size;

	hls::stream< ap_uint<5> > sizeStream1("sizeStream1");
#pragma HLS STREAM variable=sizeStream1 depth=10 dim=1
#pragma HLS RESOURCE variable=sizeStream1 core=FIFO_SRL

	hls::stream< ap_uint<5> > sizeStream2("sizeStream2");
#pragma HLS STREAM variable=sizeStream2 depth=10 dim=1
#pragma HLS RESOURCE variable=sizeStream2 core=FIFO_SRL

	hls::stream< ap_uint<5> > sizeStream3("sizeStream3");
#pragma HLS STREAM variable=sizeStream3 depth=10 dim=1
#pragma HLS RESOURCE variable=sizeStream3 core=FIFO_SRL

	ap_uint<5> size1, size2, size3;

    hls::stream<X_TYPE>  xStream("xStream");
#pragma HLS RESOURCE variable=xStream core=FIFO_SRL
    hls::stream<Y_TYPE>  yStream("yStream");
#pragma HLS RESOURCE variable=yStream core=FIFO_SRL
    hls::stream< ap_uint<1> >  polStream("polStream");
#pragma HLS RESOURCE variable=polStream core=FIFO_SRL
    hls::stream< ap_uint<TS_TYPE_BIT_WIDTH> > tsStream("tsStream");
#pragma HLS RESOURCE variable=tsStream core=FIFO_SRL
    hls::stream< sliceIdx_t > idxStream("idxStream");
#pragma HLS RESOURCE variable=idxStream core=FIFO_SRL

    hls::stream< ap_uint<96> > pktEventDataStream("pktEventDataStream");
#pragma HLS STREAM variable=pktEventDataStream depth=10 dim=1
#pragma HLS RESOURCE variable=pktEventDataStream core=FIFO_SRL

	hls::stream< ap_uint<2> >  stageOutStream("stageOutStream");
#pragma HLS STREAM variable=stageOutStream depth=10 dim=1
#pragma HLS RESOURCE variable=stageOutStream core=FIFO_SRL

	hls::stream< ap_uint<1> > stageCornerStream("stageCornerStream");

	ap_uint<5*OUTER_SIZE> idxDataWide;
    ap_uint<5> idxData[OUTER_SIZE];
    ap_uint<1> isStageCorner;

    ap_uint<2> stageOut;

//    glConfig = config;

    GetData: preProcessStream(xStreamIn, yStreamIn, polStreamIn, idxStreamIn, tsStreamIn, xStream, yStream, polStream, idxStream, tsStream, pktEventDataStream);
    Processing:
	{
//		initStageInterleaveStream(&stageOut);
//		rwSAEPerfectLoopStream<2>(xStream, yStream, tsStream, stageOutStream, outer, &size);
//		rwSAEStreamV2<2>(xStream, yStream, tsStream, stageOutStream, outer, sizeStream1, sizeStream2, sizeStream3);
		rwSAEPerfectLoopStreamV2<2>(xStream, yStream, polStream, idxStream, tsStream, inStream, &size1, &size2);
//		sizeStream1 >> size1;
//		convertInterface<4>(outer, size1, inStream);
//		sizeStream2 >> size2;
//		sortedIdxStream<5>(inStream, size1, idxData);
		sortedIdxStreamV3<4>(inStream, size1, &idxDataWide);
		getIdxDataRegion:
		{
//			for(int i = 0; i < OUTER_SIZE; i++)
//			{
//	#pragma HLS UNROLL
//				idxData[i] = idxDataWide.range(4, 0);
//				idxDataWide = idxDataWide >> 5;
//			}
//			sizeStream2 >> size2;
		}
		checkIdxGeneralV3<4>(idxDataWide, size2, &isStageCorner);   // If resource is not enough, decrease this number to increase II a little.
		feedbackInterleaveStream(isStageCorner, stageCornerStream);
	}
    Output: combineOutputStream(pktEventDataStream, stageCornerStream, xStreamOut, yStreamOut, polStreamOut, tsStreamOut, isFinalCornerStream);
//	*status = glStatus;
}

