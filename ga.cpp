/*----------------------------------------------------------
	遺伝的アルゴリズム4 実習2（時間割問題）
	作成： 2019.4.17
----------------------------------------------------------*/
/*
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "ga.h"



#define LG	20	// 遺伝子長（整数を扱う場合は32以下）
#define NP	10000   	// 個体数（偶数）
#define NG	400	// 世代数
#define NE	1000	// エリート数（偶数）
#define	RM	0.1	// 突然変異率（0～1）
#define	RC	0.7	// 交叉率（0～1）
#define	AF	2.0	// 適応度の調整（f(x)に依存）

#define MaxC 3200 //ナップサックの重量制限
#define MinSubject 10 //最小履修科目数
#define MaxSubject 12 //最大履修科目数
#define MinSpecialitySubject 5 //最小専門履修科目数
#define MinHumanSubject 2 //最小人間科学履修科目数
#define MinEnglishSubject 1 //最小英語履修科目数
#define TotalRequiredSubject 3 //合計必修履修科目数

#define NRAND(n)	(rand() % (n))	// 0 ～ n-1 の乱数（マクロ）
#define DRAND()		( ((double)rand() + 1.0) / ((double)RAND_MAX + 2.0) )
#define SWAP(type, x, y) do {type t; t=x; x=y; y=t;} while (0)	// 交換


void ga::get(std::vector<int> gene_number,std::vector<int> gene, std::ostream &out){

    //---------- 変数宣言 ----------

    char G[NP][LG];    // 個体集合
    char GW[NP][LG];    // 個体集合（作業領域）

    double Fx[NP];        // 関数値
    double RB[NP];        // ルーレットの境界
    int ord[NP];    // 並び替え用

    int i, j, k, n;    // 作業変数
    int i1, i2, j1, j2;    // 作業変数
    int ig;        // 世代変数
    int TotalSubject[20] = {};
    int x = 0;        // 作業変数


    int Satis[20][3] = {{80,0,0},{80,80,70},{70,65,0},{0,0,0},
                        {70,70,0},{70,0,0},{60,75,0},{75,0,0},
                        {90,65,0},{70,0,0},{60,0,0},{90,0,0},
                        {90,0,0},{80,50,0},{80,0,0},{80,0,0},
                        {70,0,0},{60,0,0},{0,0,0},{0,0,0}};
    int Cate[20][3] = {{1,0,0},{1,2,3},{1,3,0},{0,0,0}, // 1:専門科目 2:人間科目 3:英語 4:専門科目必修
                       {3,1,0},{4,0,0},{3,2,0},{2,0,0},
                       {3,2,0},{2,0,0},{1,0,0},{1,0,0},
                       {1,0,0},{1,3,0},{4,0,0},{4,0,0},
                       {2,0,0},{1,0,0},{0,0,0},{0,0,0}};

    srand((unsigned)time(NULL));	// 毎回異なる乱数
    int size = (int)gene_number.size();
    for(int i=0;i<size; i++){
        printf("%d",gene_number[i]);
    }


    //---------- 初期個体生成 ----------
    for (i = 0; i < NP; i++) {
        for (j = 0; j < LG; j++) {
            G[i][j] = NRAND(4);    // 0:未履修 1:1科目目　2:2科目目   3:3科目目 (x科目目に該当がなければ未履修に変換)
            printf("%d", G[i][j]);
        }
        printf("\t");
    }

//========== 世代ループ (開始) ==========
    for (ig = 0; ig <= NG; ig++) {
        int TotalSatis[NP] = {}, TotalCate[NP][5] = {}; // 1:専門科目+専門科目必修 2:人間科目 3:英語 4:専門科目必修
//---------- 関数値と適応度の計算 ----------
        for (i = 0; i < NP; i++) {

//--------------------------------------------------
//	この部分は対象とする問題に応じて作る
//--------------------------------------------------
            for (k = 0; k < LG; k++) {
                if (G[i][k] > 0) {
                    x++;
                    TotalSatis[i] += Satis[k][G[i][k]-1];
                    if(Cate[k][G[i][k]-1] == 4) TotalCate[i][4] += 1;
                    if (Cate[k][G[i][k]-1] == 1 || Cate[k][G[i][k]-1] == 4) TotalCate[i][1] += 1;
                    if (Cate[k][G[i][k]-1] == 2) TotalCate[i][2] += 1;
                    if (Cate[k][G[i][k]-1] == 3) TotalCate[i][3] += 1;
                    if (Cate[k][G[i][k]-1] == 0) G[i][k] = 0;
                }else{

                    if(k%4 == 1 && G[i][k-1] != 0 && G[i][k+1] != 0){
                        TotalSatis[i] -= 15;
                    }
                    else if(k%4 == 2 && G[i][k-1] != 0 && G[i][k+1] != 0){
                        TotalSatis[i] -= 15;
                    }
                    else if(k%4 == 2 && G[i][k-2] != 0 && G[i][k-1] == 0 && G[i][k+1] != 0){
                        TotalSatis[i] -= 30;
                    }
                }
                if(k%4 == 3 && x == 0){
                    if( (k+1)/4 == 1 || (k+1)/4 == 5){
                        TotalSatis[i] += 60;
                    }else{
                        TotalSatis[i] += 50;
                    }
                }
                if(k%4 == 3) x = 0;
            }

            //
            TotalSubject[i] = TotalCate[i][1] + TotalCate[i][2] + TotalCate[i][3];
            if( !(TotalSubject[i] >= MinSubject && TotalSubject[i] <= MaxSubject) ) TotalSatis[i] = 0;
            if( TotalCate[i][1] < MinSpecialitySubject || TotalCate[i][2] < MinHumanSubject || TotalCate[i][3] < MinEnglishSubject || TotalCate[i][4] < TotalRequiredSubject ){
                TotalSatis[i] = 0;
            }

//--------------------------------------------------
        }
        printf("\n");
        //---------- 結果表示 ----------
        Disp(ig, G[0], LG, TotalSatis[0], TotalCate[0], TotalSubject[0]);

        //===== 世代ループ (実行終了) ====
        if (ig == NG) break;

        //---------- 適応度で並び替えて作業領域にコピー (GW <- G) ---
        for (i = 0; i < NP; i++) ord[i] = i;	// 並び替え順序

        for (i = NP-1; i >= 1; i--) {		// 降順にバブルソート
            for (j = i-1; j >= 0; j--) {
                if (TotalSatis[j] < TotalSatis[i]) {
                    SWAP(int, TotalSatis[j], TotalSatis[i]);
                    SWAP(int, ord[j], ord[i]);
                }
            }
        }

        for (i = 0; i < NP; i++) {		// 並び替えてコピー
            j = ord[i];
            for (k = 0; k < LG; k++) GW[i][k] = G[j][k];
        }


        //---------- 次世代個体生成 (G <- GW) ----------

        //----- ルーレット境界生成 -----
        RB[0] = TotalSatis[0];
        for (i = 0; i < NP; i++) RB[i] = RB[i-1] + TotalSatis[i];
        for (i = 0; i < NP; i++) RB[i] /= RB[NP-1];

        //----- エリート保存 -----
        for (i = 0; i < NE; i++) {
            for (k = 0; k < LG; k++) {
                G[i][k] = GW[i][k];
            }
        }

        //----- エリート以外の個体生成 -----
        for (n = 0; n < (NP - NE)/2; n++) {

            //--- 選択 ---
            i1 = i2 = RouSel(RB);
            while (i1 == i2) i2 = RouSel(RB);

            //--- 交叉 ---
            j1 = 2 * (n + 1);
            j2 = 2 * (n + 1) + 1;
            if (DRAND() < RC) {

                //--- 一様交叉 ---//
                for (k = 0; k < LG; k++)
                    if (NRAND(2)) SWAP(char, G[j1][k],G[j2][k]);
            }
        }

        //----- 突然変異 -----
        for (i = NE + 1; i < NP; i++) {
            for (k = 0; k < LG; k++) {
                if (DRAND() < RM) G[i][k] = NRAND(4);
            }
        }
    }
    //========== 世代ループ (記述終了) ==========

}


//---------- 個体の選択 ----------

int	RouSel(double y[]) {
    int i;
    double r;

    i = 0;
    r = DRAND();
    while(r > y[i]) i++;

    return i;
}

//---------- 変数 → 画面 ----------

int Disp(int i, char x[], int L, int Satis, int Cate[], int TotalSubject) {
    int k;

    printf("%d,\t\"", i);
    for (k = 0; k < L; k++){
        printf("%d", x[k]);
        if(k%4 == 3) printf(" ");
    }
    printf("\",\t%d,\t%d,\t%d,\t%d,\t%d,\t%d\n", Satis,Cate[1],Cate[2],Cate[3],Cate[4],TotalSubject);
    return 0;
}*/