#include <iostream>
#include <vector>
#include <complex>
#include <algorithm>
#include <iomanip>
#include <cmath>

#define MINIMP3_IMPLEMENTATION

#include "minimp3.h"
#include "minimp3_ex.h"

using namespace std;

const int BLOCK_SIZE = 4096;
const int STEP_SIZE = 1024;

void setUpIoOptimizations();

vector<short> decodeMp3File(const char *fileName);

vector<double> getMaximumValueByBlock(const vector<short> &decodedMp3);

void printResult(vector<double> &maxValuesInBlocks);

void applyFFT(vector<complex<double>> &a);

vector<complex<double>> applyHannWindow(const vector<short> &block);

vector<short> prepareBlockForHannWindow(const vector<short> &decodedMp3, unsigned int l, unsigned int r);

complex<double> getMaxElementInBlock(vector<complex<double>> &blockAfterHannWindowApplying);

double
getMaxAbsValueInBlock(const vector<short> &decodedMp3, unsigned int leftBoundOfBlock, unsigned int rightBoundOfBlock);

int main() {
    setUpIoOptimizations();

    vector<short> decodedMp3 = decodeMp3File("input.mp3");

    vector<double> maxValuesInBlocks = getMaximumValueByBlock(decodedMp3);

    printResult(maxValuesInBlocks);
    return 0;
}

void printResult(vector<double> &maxValuesInBlocks) {
    for (double cnt: maxValuesInBlocks) {
        cout << fixed << setprecision(20) << cnt << "\n";
    }
}

vector<double> getMaximumValueByBlock(const vector<short> &decodedMp3) {
    vector<double> maxValuesInBlocks;

    unsigned leftBoundOfBlock = 0;
    unsigned rightBoundOfBlock = BLOCK_SIZE;

    while (rightBoundOfBlock <= decodedMp3.size()) {
        double maxAbsValueInBlock = getMaxAbsValueInBlock(decodedMp3, leftBoundOfBlock, rightBoundOfBlock);

        maxValuesInBlocks.push_back(maxAbsValueInBlock);

        leftBoundOfBlock += STEP_SIZE;
        rightBoundOfBlock += STEP_SIZE;
    }

    return maxValuesInBlocks;
}

double
getMaxAbsValueInBlock(const vector<short> &decodedMp3, unsigned int leftBoundOfBlock, unsigned int rightBoundOfBlock) {
    vector<short> blockForHannWindow = prepareBlockForHannWindow(decodedMp3, leftBoundOfBlock, rightBoundOfBlock);
    vector<complex<double>> blockAfterHannWindowApplying = applyHannWindow(blockForHannWindow);
    applyFFT(blockAfterHannWindowApplying);

    complex<double> maxElementInBlock = getMaxElementInBlock(blockAfterHannWindowApplying);

    return abs(maxElementInBlock);
}

complex<double> getMaxElementInBlock(vector<complex<double>> &blockAfterHannWindowApplying) {
    return *max_element(
            blockAfterHannWindowApplying.begin(),
            blockAfterHannWindowApplying.end(),
            [](complex<double> a, complex<double> b) {
                return abs(a) < abs(b);
            });
}

vector<short> prepareBlockForHannWindow(const vector<short> &decodedMp3, unsigned int l, unsigned int r) {
    vector<short> blockForHannWindow;
    for (unsigned i = l; i < r; ++i) {
        blockForHannWindow.push_back(decodedMp3[i]);
    }
    return blockForHannWindow;
}


void setUpIoOptimizations() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    cout.tie(nullptr);
}

vector<short> decodeMp3File(const char *fileName) {
    mp3dec_t mp3decoder;
    mp3dec_file_info_t decodedFileInfo;
    mp3dec_load(&mp3decoder, fileName, &decodedFileInfo, nullptr, nullptr);
    vector<short> decodedFile(decodedFileInfo.buffer, decodedFileInfo.buffer + decodedFileInfo.samples);
    return decodedFile;
}

int reverse(int num, int lgN) {
    int res = 0;
    for (int i = 0; i < lgN; i++) {
        if (num & (1 << i))
            res |= 1 << (lgN - 1 - i);
    }
    return res;
}

vector<complex<double>> applyHannWindow(const vector<short> &block) {
    vector<double> calculatedValues;
    for (unsigned elem = 0; elem < BLOCK_SIZE; ++elem) {
        double multiplier = 0.5 * (1 - cos(2 * M_PI * elem / (BLOCK_SIZE - 1)));
        calculatedValues.push_back(block[elem] * multiplier);
    }
    vector<complex<double>> blockAfterHannWindow(calculatedValues.begin(), calculatedValues.end());
    return std::move(blockAfterHannWindow);
}

void applyFFT(vector<complex<double>> &a) {
    int n = (int) a.size();

    int lgN = 0;
    while ((1 << lgN) < n)
        lgN++;

    for (int i = 0; i < n; i++) {
        if (i < reverse(i, lgN))
            swap(a[i], a[reverse(i, lgN)]);
    }

    for (int lengthOfStep = 2; lengthOfStep <= n; lengthOfStep <<= 1) {
        double mainRoot = 2 * M_PI / lengthOfStep;
        complex<double> nextRootMultiplier(cos(mainRoot), sin(mainRoot));
        for (int i = 0; i < n; i += lengthOfStep) {
            complex<double> currentRoot(1);
            for (int j = 0; j < lengthOfStep / 2; j++) {
                complex<double> u = a[i + j];
                complex<double> v = a[i + j + lengthOfStep / 2] * currentRoot;
                a[i + j] = u + v;
                a[i + j + lengthOfStep / 2] = u - v;
                currentRoot *= nextRootMultiplier;
            }
        }
    }
}
