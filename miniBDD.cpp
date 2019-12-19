#include <cassert>

#include "miniBDD.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <fstream>
using std::endl;
using std::ofstream;


//#define LG	20	// 遺伝子長（整数を扱う場合は32以下）
#define NP	2000   	// 個体数（偶数）
#define NG	300	// 世代数
#define NE	300	// エリート数（偶数）
#define	RM	0.8	// 突然変異率（0～1）
#define	RC	0.7	// 交叉率（0～1）
#define	AF	10.0	// 適応度の調整（f(x)に依存）
#define ER  0.8  //許容誤差

#define NRAND(n)	(rand() % (n))	// 0 ～ n-1 の乱数（マクロ）
#define DRAND()		( ((double)rand() + 1.0) / ((double)RAND_MAX + 2.0) )
#define SWAP(type, x, y) do {type t; t=x; x=y; y=t;} while (0)	// 交換



void BDD::clear()
{
  if(node!=NULL)
  {
    node->remove_reference();
    node=NULL;
  }
}

void BDDnode::remove_reference()
{
  assert(reference_counter!=0);
  
  reference_counter--;

  if(reference_counter==0 && node_number>=2)
  {
    miniBDD_mgr::reverse_keyt reverse_key(var, low, high);
    mgr->reverse_map.erase(reverse_key);
    low.clear();
    high.clear();
  } 
}

BDD miniBDD_mgr::Var(const std::string &label)
{
  var_table.push_back(var_table_entryt(label));
  true_bdd.node->var=var_table.size()+1;
  false_bdd.node->var=var_table.size()+1;
  return mk(var_table.size(), false_bdd, true_bdd);
}

void miniBDD_mgr::DumpDot(std::ostream &out, bool suppress_zero) const
{
  out << "digraph BDD {" << std::endl;

  out << "center = true;" << std::endl;
  
  out << "{ rank = same; { node [style=invis]; \"T\" };" << std::endl;
  
  if(!suppress_zero)
    out << "  { node [shape=box,fontsize=24]; \"0\"; }" << std::endl;
    
  out << "  { node [shape=box,fontsize=24]; \"1\"; }" << std::endl
      << "}" << std::endl
      << std::endl;
      
  for(unsigned v=0; v<var_table.size(); v++)
  {
    out << "{ rank=same; "
           "{ node [shape=plaintext,fontname=\"Times Italic\",fontsize=24] \" "
        << var_table[v].label
        << " \" }; ";

    forall_nodes(u)
      if(u->var==(v+1) && u->reference_counter!=0)
        out << '"' << u->node_number << "\"; ";
    
    out << "}" << std::endl;
  }

  out << std::endl;

  out << "{ edge [style = invis];";

  for(unsigned v=0; v<var_table.size(); v++)
    out << " \" " << var_table[v].label
        << " \" ->";
  
  out << " \"T\"; }" << std::endl;
  
  out << std::endl;

  forall_nodes(u)
  {
    if(u->reference_counter==0) continue;
    if(u->node_number<=1) continue;

    if(!suppress_zero || u->high.node_number()!=0)
      out << '"' << u->node_number << '"' << " -> "
          << '"' << u->high.node_number() << '"'
          << " [style=solid,arrowsize=\".75\"];"
          << std::endl;
        
    if(!suppress_zero || u->low.node_number()!=0)
      out << '"' << u->node_number << '"' << " -> "
          << '"' << u->low.node_number() << '"'
          << " [style=dashed,arrowsize=\".75\"];"
          << std::endl;

    out << std::endl;
  }
  
  out << "}" << std::endl;
}

void miniBDD_mgr::DumpTikZ(
  std::ostream &out,
  bool suppress_zero,
  bool node_numbers) const
{
  out << "\\begin{tikzpicture}[node distance=1cm]" << std::endl;
  
  out << "  \\tikzstyle{BDDnode}=[circle,draw=black,"
         "inner sep=0pt,minimum size=5mm]" << std::endl;

  for(unsigned v=0; v<var_table.size(); v++)
  {
    out << "  \\node[";

    if(v!=0)
      out << "below of=v" << var_table[v-1].label;

    out << "] (v" << var_table[v].label << ") {$\\mathit{"
        << var_table[v].label << "}$};" << std::endl;

    unsigned previous=0;

    forall_nodes(u)
    {
      if(u->var==(v+1) && u->reference_counter!=0)
      {
        out << "  \\node[xshift=0cm, BDDnode, ";

        if(previous==0)
          out << "right of=v" << var_table[v].label;
        else
          out << "right of=n" << previous;

        out << "] (n" << u->node_number << ") {";
        if(node_numbers) out << "\\small $" << u->node_number << "$";
        out << "};" << std::endl;
        previous=u->node_number;
      }
    }
    
    out << std::endl;
  }

  out << std::endl;

  out << "  % terminals" << std::endl;
  out << "  \\node[draw=black, style=rectangle, below of=v"
      << var_table.back().label
      << ", xshift=1cm] (n1) {$1$};" << std::endl;
    
  if(!suppress_zero)
    out << "  \\node[draw=black, style=rectangle, left of=n1] (n0) {$0$};" << std::endl;

  out << std::endl;

  out << "  % edges" << std::endl;
  out << std::endl;

  forall_nodes(u)
  {
    if(u->reference_counter!=0 && u->node_number>=2)
    {
      if(!suppress_zero || u->low.node_number()!=0)
        out << "  \\draw[->,dashed] (n" << u->node_number << ") -> (n"
            << u->low.node_number() << ");" << std::endl;
          
      if(!suppress_zero || u->high.node_number()!=0)
        out << "  \\draw[->] (n" << u->node_number << ") -> (n"
            << u->high.node_number() << ");" << std::endl;
    }
  }

  out << std::endl;
  
  out << "\\end{tikzpicture}" << std::endl;
}

bool equal_fkt(bool x, bool y)
{
  return x==y;
}

BDD BDD::operator ==(const BDD &other) const
{
  return apply(equal_fkt, *this, other);
}

bool xor_fkt(bool x, bool y)
{
  return x ^ y;
}

BDD BDD::operator ^(const BDD &other) const
{
  return apply(xor_fkt, *this, other);
}

BDD BDD::operator !() const
{
  return node->mgr->True() ^ *this;
}

bool and_fkt(bool x, bool y)
{
  return x && y;
}

BDD BDD::operator &(const BDD &other) const
{
  return apply(and_fkt, *this, other);
}

bool or_fkt(bool x, bool y)
{
  return x || y;
}

BDD BDD::operator |(const BDD &other) const
{
  return apply(or_fkt, *this, other);
}

miniBDD_mgr::miniBDD_mgr()
{
  // add true/false nodes
  nodes.push_back(BDDnode(this, 0, 0, BDD(), BDD()));
  false_bdd=BDD(&nodes.back());
  nodes.push_back(BDDnode(this, 1, 1, BDD(), BDD()));
  true_bdd=BDD(&nodes.back());
}

miniBDD_mgr::~miniBDD_mgr()
{
}

BDD miniBDD_mgr::mk(unsigned var, const BDD &low, const BDD &high)
{
  assert(var<=var_table.size());

  if(low.node_number()==high.node_number())
    return low;
  else
  {
    reverse_keyt reverse_key(var, low, high);
    reverse_mapt::const_iterator it=reverse_map.find(reverse_key);

    if(it!=reverse_map.end())
      return BDD(it->second);
    else
    {
      unsigned new_number=nodes.back().node_number+1;
      nodes.push_back(BDDnode(this, var, new_number, low, high));
      reverse_map[reverse_key]=&nodes.back();
      return BDD(&nodes.back());
    }
  }
}

bool operator < (const miniBDD_mgr::reverse_keyt &x,
                 const miniBDD_mgr::reverse_keyt &y)
{
  if(x.var<y.var) return true;
  if(x.var>y.var) return false;
  if(x.low<y.low) return true;
  if(x.low>y.low) return false;
  return x.high<y.high;
}

void miniBDD_mgr::DumpTable(std::ostream &out) const
{
  out << "\\# & \\mathit{var} & \\mathit{low} &"
         " \\mathit{high} \\\\\\hline" << std::endl;

  forall_nodes(it)
  {
    out << it->node_number << " & ";

    if(it->node_number==0 || it->node_number==1)
      out << it->var << " & & \\\\";
    else if(it->reference_counter==0)
      out << "- & - & - \\\\";
    else
      out << it->var << "\\," << var_table[it->var-1].label << " & "
          << it->low.node_number() << " & " << it->high.node_number()
          << " \\\\";
          
    if(it->node_number==1) out << "\\hline";

    out << " % " << it->reference_counter << std::endl;
  }
}

void miniBDD_mgr::DumpNumber(std::ostream &out) const
{
  std::vector<int> gene;
  std::vector<int> gene_number;
  std::vector<int> gene_var;
  std::vector<int> gene_high;
  std::vector<int> gene_low;
  std::vector<int> jump_gene;


    int LG = 0, LastNode = -1;
    char OriginalFlug = 0;

    forall_nodes(it)
    {
        LastNode++;

        if(!(it->node_number==0 || it->node_number==1) && !(it->reference_counter==0)){
            LG++;
            gene.push_back(1);
            gene_number.push_back(it->node_number);
            gene_var.push_back(it->var);
            gene_high.push_back(it->high.node_number());
            gene_low.push_back(it->low.node_number());

        }
    }

    char G[NP][LG];    // 個体集合
    char GW[NP][LG];    // 個体集合（作業領域）

    double Fx[NP];        // 関数値
    double RB[NP];        // ルーレットの境界
    int ord[NP];    // 並び替え用

    int i, j, k, n,x;    // 作業変数
    int i1, i2, j1, j2;    // 作業変数
    int ig;        // 世代変数
    int OriginalTrue = 0;
    int OriginalLow[LG], OriginalHigh[LG];
    int Current_gene = 0, Current_LH = 0, Current_count = 0;

    srand((unsigned)time(NULL));	// 毎回異なる乱数
    //---------- 初期個体生成 ----------
    for (i = 0; i < NP; i++) {
        for (j = 0; j < LG; j++) {
            G[i][j] = NRAND(2);    // 0:未履修 1:1科目目　2:2科目目   3:3科目目 (x科目目に該当がなければ未履修に変換)
            printf("%d", G[i][j]);
        }
        printf("\t");
    }
    printf("\nNumber\tGene\t  True  ErrorRate\tFitness\t\tAbEr");
//========== 世代ループ (開始) ==========
    for (ig = 0; ig <= NG; ig++) {
        int TotalTrue[NP] = {}, TotalGene[NP] = {}; // 生根の数
        int ScoreLow[NP][LG], ScoreHigh[NP][LG];
        for (i = 0; i < NP; i++) {
            for (j = 0; j < LG; j++) {
                ScoreLow[i][j] = 0; ScoreHigh[i][j] = 0;
            }
        }
        double Fitness[NP] = {}, AbEr[NP] = {}, ErrorRate[NP] = {};

//---------- 関数値と適応度の計算 ----------
        for (i = 0; i < NP; i++) {
            int k = 0; int n = 0, varMax = 0, weight = 0;
            int cloneFlag = 0;
//--------------------------------------------------
//	この部分は対象とする問題に応じて作る
//--------------------------------------------------
            forall_nodes(it) {
                if(it->node_number==0 || it->node_number==1){varMax = it->var;}
                else if(it->reference_counter==0){}
                else {
                    weight = varMax - (it->var) - 1;
                    if(OriginalFlug < LG){
                        OriginalLow[k] = 0;
                        OriginalHigh[k] = 0;
                        if(it->low.node_number() == 1){
                            //OriginalTrue += pow(2,weight);
                            OriginalLow[k] += pow(2,weight);
                        }
                        else if(it->high.node_number() == 1){
                            //OriginalTrue += pow(2,weight);
                            OriginalHigh[k] += pow(2,weight);
                        }
                        if(it->low.node_number() != 0 && it->low.node_number() != 1){
                            x =it->low.node_number();
                            for(n=0;gene_number[n]<=x;n++) {
                                if(gene_number[n] == x){
                                    OriginalLow[k] += OriginalLow[n] + OriginalHigh[n];
                                }
                            }
                        }
                        if(it->high.node_number() != 0 && it->high.node_number() != 1){
                            x =it->high.node_number();
                            for(n=0;gene_number[n]<=x;n++) {
                                if(gene_number[n] == x){
                                    OriginalHigh[k] += OriginalLow[n] + OriginalHigh[n];
                                }
                            }
                        }
                        if(LastNode == it->node_number){OriginalTrue += OriginalLow[k] + OriginalHigh[k];}
                        OriginalFlug++;
                    }
                    if(G[i][k]==1) {
                        if(it->low.node_number() == 1){
                            ScoreLow[i][k] += pow(2,weight);
                            TotalTrue[i] += pow(2,weight);
                        }else if(!(it->low.node_number() == 0)){
                            x =it->low.node_number();
                            for(n=0;gene_number[n] <= x;n++) {
                                if(gene_number[n] == x && G[i][n] == 0){
                                    ScoreLow[i][k] += pow(2, weight);
                                    TotalTrue[i] += pow(2,weight);
                                }
                            }
                        }
                        if(it->high.node_number() == 1){
                            ScoreHigh[i][k] += pow(2,weight);
                            TotalTrue[i] += pow(2,weight);
                        }else if(!(it->high.node_number() == 0)){
                            x =it->high.node_number();
                            for(n=0;gene_number[n] <= x;n++) {
                                if(gene_number[n] == x && G[i][n] == 0){
                                    ScoreHigh[i][k] += pow(2, weight);
                                    TotalTrue[i] += pow(2,weight);
                                }
                            }
                        }

                        //if(it->high.node_number() == 1)Fitness[i] += 2^k;
                        TotalGene[i]++;
                    }else{
                        ScoreLow[i][k] *= 2;
                        ScoreHigh[i][k] *= 2;
                        TotalTrue[i] *=2;
                    }
                    k++;
                }
            }
            for(k=0;k<LG;k++){
                if(OriginalLow[k] != ScoreLow[i][k] && G[i][k] == 1){
                    ErrorRate[i] += abs(OriginalLow[k] - ScoreLow[i][k]);
                }

                if(OriginalHigh[k] != ScoreHigh[i][k] && G[i][k] == 1){
                    ErrorRate[i] += abs(OriginalHigh[k] - ScoreHigh[i][k]);
                }
            }
                if(abs(OriginalTrue-TotalTrue[i]) == 0 || TotalGene[i] == 0){
                Fitness[i] = 0;
                ErrorRate[i] = 0;
            }else{
                ErrorRate[i] += abs(double(TotalTrue[i]) - double(OriginalTrue)) / OriginalTrue / AF;
                Fitness[i] = AF / ErrorRate[i] / TotalGene[i];
                AbEr[i] = abs(double(TotalTrue[i]) - double(OriginalTrue));
                if(AbEr[i]/pow(2,varMax-1) <= ER){
                    Fitness[i] += 10;
                }else{

                }
            }

//--------------------------------------------------
        }

        printf("\n");

        //---------- 結果表示 ----------
        //Disp(ig, G[0], LG, TotalSatis[0], TotalCate[0], TotalSubject[0]);
        printf("%d,\t", ig);
        for (k = 0; k < LG; k++){
            printf("%d", G[0][k]);
        }
        printf("    %d\t%lf\t%lf\t%lf",TotalTrue[0],ErrorRate[0],Fitness[0],AbEr[0]);


        //---------- 適応度で並び替えて作業領域にコピー (GW <- G) ---
        for (i = 0; i < NP; i++) ord[i] = i;	// 並び替え順序

        for (i = NP-1; i >= 1; i--) {		// 降順にバブルソート
            for (j = i-1; j >= 0; j--) {
                if (Fitness[j] < Fitness[i]) {
                    SWAP(double , Fitness[j], Fitness[i]);
                    SWAP(int, ord[j], ord[i]);
                    SWAP(double,AbEr[j],AbEr[i]);
                    SWAP(double,TotalTrue[j],TotalTrue[i]);
                    SWAP(double,ErrorRate[j],ErrorRate[i]);
                }
            }
        }


        //===== 世代ループ (実行終了) ====
        if (ig == NG) break;

        for (i = 0; i < NP; i++) {		// 並び替えてコピー
            j = ord[i];
            for (k = 0; k < LG; k++) GW[i][k] = G[j][k];
        }


        //---------- 次世代個体生成 (G <- GW) ----------

        //----- ルーレット境界生成 -----
        RB[0] = Fitness[0];
        for (i = 0; i < NP; i++) RB[i] = RB[i-1] + Fitness[i];
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
                if (DRAND() < RM) G[i][k] = NRAND(2);
            }
        }
    }
    //========== 世代ループ (記述終了) ==========

    ofstream ofs("test.csv");  // ファイルパスを指定する
    ofs << 0 << ", "<< 2 << ", " << endl;
    ofs << 3 << ", "<< 4 << ", " << endl;
    ofs << 5 << ", "<< 6 << ", " << endl;
    Current_gene = LG-1;
    printf("\n");
    if(gene[Current_gene]==1){
        //ALL_LOW------------------------------------------------------------
        if(gene_low[Current_gene]==0){
            for(k=0; k<pow(2,gene_var[0]-gene_var[Current_gene]); k++){
                Current_count++;
                printf("0");
            }
        }else if(gene_low[Current_gene] == 1){
            for(k=0; k<pow(2,gene_var[0]-gene_var[Current_gene]); k++){
                Current_count++;
                printf("1");
            }
        }else{
            for(n=0;gene_number[n]<gene_low[Current_gene];n++) {
            }
            Current_gene = n;
            while(gene_low[Current_gene]+gene_high[Current_gene] > 2){
                Current_gene = n;
                //Low
                if(gene_low[Current_gene]==0){
                    for(k=0; k<pow(2,gene_var[0]-gene_var[Current_gene]); k++){
                        Current_count++;
                        printf("0");
                    }
                }else if(gene_low[Current_gene] == 1){
                    for(k=0; k<pow(2,gene_var[0]-gene_var[Current_gene]); k++){
                        Current_count++;
                        printf("1");
                    }
                }else{
                    for(n=0;gene_number[n]<gene_low[Current_gene];n++) {
                    }
                    jump_gene.push_back(Current_gene);
                    continue;
                }

                //High
                if(gene_high[Current_gene]==0){
                    for(k=0; k<pow(2,gene_var[0]-gene_var[Current_gene]); k++){
                        Current_count++;
                        printf("0");
                    }
                }else if(gene_high[Current_gene] == 1){
                    for(k=0; k<pow(2,gene_var[0]-gene_var[Current_gene]); k++){
                        Current_count++;
                        printf("1");
                    }
                }else{
                    for(n=0;gene_number[n]<gene_high[Current_gene];n++) {
                    }
                }

                //飛ばした根の評価
                x = Current_gene;
                while(!jump_gene.empty()){
                    Current_gene = jump_gene.back();
                    jump_gene.pop_back();
                    if(gene_high[Current_gene]==0){
                        for(k=0; k<pow(2,gene_var[0]-gene_var[Current_gene]); k++){
                            Current_count++;
                            printf("0");
                        }
                    }else if(gene_high[Current_gene] == 1){
                        for(k=0; k<pow(2,gene_var[0]-gene_var[Current_gene]); k++){
                            Current_count++;
                            printf("1");
                        }
                    }else{
                        for(n=0;gene_number[n]<gene_high[Current_gene];n++) {
                        }
                        x = Current_gene;
                    }
                    Current_gene = x;
                }
            }
        }

        //ALL_HIGH------------------------------------------------------------
        Current_gene = LG-1;
        if(gene_high[Current_gene]==0){
            for(k=0; k<pow(2,gene_var[0]-gene_var[Current_gene]); k++){
                Current_count++;
                printf("0");
            }
        }else if(gene_high[Current_gene] == 1){
            for(k=0; k<pow(2,gene_var[0]-gene_var[Current_gene]); k++){
                Current_count++;
                printf("1");
            }
        }else{
            for(n=0;gene_number[n]<gene_high[Current_gene];n++) {
            }
            Current_gene = n;
            while(gene_low[Current_gene]+gene_high[Current_gene] > 2){
                Current_gene = n;
                //Low
                if(gene_low[Current_gene]==0){
                    for(k=0; k<pow(2,gene_var[0]-gene_var[Current_gene]); k++){
                        Current_count++;
                        printf("0");
                    }
                }else if(gene_low[Current_gene] == 1){
                    for(k=0; k<pow(2,gene_var[0]-gene_var[Current_gene]); k++){
                        Current_count++;
                        printf("1");
                    }
                }else{
                    for(n=0;gene_number[n]<gene_low[Current_gene];n++) {
                    }
                    jump_gene.push_back(Current_gene);
                    continue;
                }

                //High
                if(gene_high[Current_gene]==0){
                    for(k=0; k<pow(2,gene_var[0]-gene_var[Current_gene]); k++){
                        Current_count++;
                        printf("0");
                    }
                }else if(gene_high[Current_gene] == 1){
                    for(k=0; k<pow(2,gene_var[0]-gene_var[Current_gene]); k++){
                        Current_count++;
                        printf("1");
                    }
                }else{
                    for(n=0;gene_number[n]<gene_high[Current_gene];n++) {
                    }
                }

                //飛ばした根の評価
                x = Current_gene;
                while(!jump_gene.empty()){
                    Current_gene = jump_gene.back();
                    jump_gene.pop_back();
                    if(gene_high[Current_gene]==0){
                        for(k=0; k<pow(2,gene_var[0]-gene_var[Current_gene]); k++){
                            Current_count++;
                            printf("0");
                        }
                    }else if(gene_high[Current_gene] == 1){
                        for(k=0; k<pow(2,gene_var[0]-gene_var[Current_gene]); k++){
                            Current_count++;
                            printf("1");
                        }
                    }else{
                        for(n=0;gene_number[n]<gene_high[Current_gene];n++) {
                        }
                        x = Current_gene;
                    }
                    Current_gene = x;
                }
            }
        }
    }





}

BDD apply(bool (*fkt)(bool x, bool y), const BDD &x, const BDD &y)
{
  assert(x.node!=NULL && y.node!=NULL);
  assert(x.node->mgr==y.node->mgr);

  miniBDD_mgr *mgr=x.node->mgr;
  
  BDD u;

  if(x.is_constant() && y.is_constant())
    u=BDD(fkt(x.is_true(), y.is_true())?mgr->true_bdd:mgr->false_bdd);
  else if(x.var()==y.var())
    u=mgr->mk(x.var(),
              apply(fkt, x.low(), y.low()),
              apply(fkt, x.high(), y.high()));
  else if(x.var()<y.var())
    u=mgr->mk(x.var(),
              apply(fkt, x.low(), y),
              apply(fkt, x.high(), y));
  else /* x.var() > y.var() */
    u=mgr->mk(y.var(),
              apply(fkt, x, y.low()),
              apply(fkt, x, y.high()));
    
  return u;
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
