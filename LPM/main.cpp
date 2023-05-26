#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string.h>

// Function to implement the KMP algorithm
template <typename T>
bool KMP(T textbeg, T textend, T patbeg, T patend, std::string &pattern, bool reverse, int &lp, std::vector<int> &lpi, int threshold)
{
    std::vector<int> next(patend - patbeg, 0);

    for (auto it = patbeg + 1; it < patend; it++)
    {
        int j = next[it - patbeg - 1];

        while (j > 0 && *(patbeg + j) != *it)
            j = next[j - 1];

        if (j > 0 || *(patbeg + j) == *it)
            next[it - patbeg] = j + 1;
    }

    bool lpiclr=false;
    int j = 0;
    for (auto it = textbeg; it < textend; it++)
    {
        if (*it == *(patbeg + j))
        {
            ++j;
            if (j>=threshold && j >= lp)
            {
                if (j > lp)
                {
                    lp = j;
                    lpi.clear();
                    lpiclr=true;
                }
                if (reverse)
                    lpi.push_back(textend - it - 1);
                else
                    lpi.push_back(it - textbeg - j + 1);
            }
        }
        else if (j > 0)
        {
            j = next[j - 1];
            it--; // since `i` will be incremented in the next iteration
        }
    }

    return lpiclr;
}

std::string get_RC(std::string &text)
{
    std::string textRC(text.size(), 'N');
    for (int i=0; i<text.size(); ++i)
    {
        int j=text.size()-i-1;
        if (text[i]=='A')
            textRC[j]='T';
        else if (text[i]=='T')
            textRC[j]='A';
        else if (text[i]=='G')
            textRC[j]='C';
        else if (text[i]=='C')
            textRC[j]='G';
    }
    return textRC;
}

// Program to implement the KMP algorithm in C++
int main(int argc, char **argv)
{
    bool reverse = (argc<2 || strcmp(argv[2], "reverse")) ? false : true;
    int threshold = argc<3 ? 20 : std::stoi(argv[3]);


    std::ifstream fin(argv[1]);
    std::string text, pattern;
    int offset;
    fin >> text;
    std::transform(text.begin(), text.end(), text.begin(), ::toupper);
    fin.close();
    std::string textRC = get_RC(text);
    
    while (std::cin >> offset)
    {
        std::cin >> pattern;
        int lp=0, Wlen;
        std::vector<int> lpi;
        if (reverse)
        {
            KMP(text.rbegin(), text.rend(), pattern.rbegin() + offset, pattern.rend(), pattern, reverse, lp, lpi, threshold);
            Wlen=lpi.size();
            if (KMP(textRC.rbegin(), textRC.rend(), pattern.rbegin() + offset, pattern.rend(), pattern, reverse, lp, lpi, threshold))
                Wlen=0;
        }
        else
        {
            KMP(text.begin(), text.end(), pattern.begin() + offset, pattern.end(), pattern, reverse, lp, lpi, threshold);
            Wlen=lpi.size();
            if (KMP(textRC.begin(), textRC.end(), pattern.begin() + offset, pattern.end(), pattern, reverse, lp, lpi, threshold))
                Wlen=0;
        }
        if (lp >= threshold)
        {
            std::cout << lp;
            for (int i = 0; i < Wlen; ++i)
                std::cout << '\t' << lpi[i] << "\t+";
            for (int i = Wlen; i < lpi.size(); ++i)
                std::cout << '\t' << text.size()-lpi[i]-lp << "\t-";
            std::cout << '\t' << pattern << '\n';
        }
    }

    return 0;
}

