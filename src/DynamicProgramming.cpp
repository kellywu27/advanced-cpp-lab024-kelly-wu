#include "project/DynamicProgramming.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>

std::string buildParenthesization(
    const std::vector<std::vector<int>>& split,
    int i,
    int j) {

    // base case
    if (i == j) {
        return "A" + std::to_string(i);
    }

    // rec case: split at k
    const int k = split[i][j];
    return "(" + buildParenthesization(split, i, k) + " x " +
           buildParenthesization(split, k + 1, j) + ")";
}

MatrixChainResult matrixChainOrder(
    const std::vector<int>& dimensions) {

    if (dimensions.size() < 2) {
        throw std::invalid_argument("At least one matrix is required");
    }
    for (int d : dimensions) {
        if (d <= 0) {
            throw std::invalid_argument("All dimensions must be positive");
        }
    }

    // n matrices; A_i is dimensions[i-1] x dimensions[i], indexed 1..n.
    const int n = static_cast<int>(dimensions.size()) - 1;

    // cost[i][j] = minimum scalar multiplications 
    std::vector<std::vector<long long>> cost(
        n + 1, std::vector<long long>(n + 1, 0));
    std::vector<std::vector<int>> split(
        n + 1, std::vector<int>(n + 1, 0));

    // evaluate chains of inc length
    for (int length = 2; length <= n; ++length) {
        for (int i = 1; i + length - 1 <= n; ++i) {
            const int j = i + length - 1;
            cost[i][j] = std::numeric_limits<long long>::max();

            for (int k = i; k < j; ++k) {
                const long long candidate =
                    cost[i][k] + cost[k + 1][j] +
                    static_cast<long long>(dimensions[i - 1]) *
                        dimensions[k] * dimensions[j];

                if (candidate < cost[i][j]) {
                    cost[i][j] = candidate;
                    split[i][j] = k;
                }
            }
        }
    }

    MatrixChainResult result;
    result.minimumCost = cost[1][n];
    result.split = split;
    result.parenthesization = buildParenthesization(split, 1, n);
    return result;
}

AlignmentResult globalSequenceAlignment(
    const std::string& first,
    const std::string& second,
    int matchScore,
    int mismatchScore,
    int gapScore) {

    const int n = static_cast<int>(first.size());
    const int m = static_cast<int>(second.size());

    // dp[i][j] = best score aligning first[0..i) with second[0..j).
    std::vector<std::vector<int>> dp(n + 1, std::vector<int>(m + 1, 0));

    for (int i = 1; i <= n; ++i) {
        dp[i][0] = i * gapScore;
    }
    for (int j = 1; j <= m; ++j) {
        dp[0][j] = j * gapScore;
    }

    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            const int diagonal = dp[i - 1][j - 1] +
                (first[i - 1] == second[j - 1] ? matchScore : mismatchScore);
            const int up = dp[i - 1][j] + gapScore;    // gap in second
            const int left = dp[i][j - 1] + gapScore;  // gap in first
            dp[i][j] = std::max({diagonal, up, left});
        }
    }

    // backtrack from the final cell, pref diagonal, then up, then left
    std::string alignedFirst;
    std::string alignedSecond;
    int i = n;
    int j = m;

    while (i > 0 || j > 0) {
        if (i > 0 && j > 0 &&
            dp[i][j] == dp[i - 1][j - 1] +
                (first[i - 1] == second[j - 1] ? matchScore : mismatchScore)) {
            alignedFirst += first[i - 1];
            alignedSecond += second[j - 1];
            --i;
            --j;
        } else if (i > 0 && dp[i][j] == dp[i - 1][j] + gapScore) {
            alignedFirst += first[i - 1];
            alignedSecond += '-';
            --i;
        } else {
            alignedFirst += '-';
            alignedSecond += second[j - 1];
            --j;
        }
    }

    // chars were appended while walking backward; reverse them
    std::reverse(alignedFirst.begin(), alignedFirst.end());
    std::reverse(alignedSecond.begin(), alignedSecond.end());

    return {dp[n][m], alignedFirst, alignedSecond};
}