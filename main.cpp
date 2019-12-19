
#include <iostream>

#include "miniBDD.h"
#include "ga.h"


void test1()
{
    miniBDD_mgr mgr;

    BDD x=mgr.Var("x");
    BDD y=mgr.Var("y");
    BDD z=mgr.Var("z");
    BDD f=(x&y&z)|(!x&!y&z);
    y.clear();
    x.clear();
    z.clear();

    //mgr.DumpDot(std::cout);
    mgr.DumpTikZ(std::cout);
}

void test2()
{
    miniBDD_mgr mgr;

    BDD a=mgr.Var("a");
    BDD b=mgr.Var("b");
    BDD c=mgr.Var("c");
    BDD d=mgr.Var("d");
    BDD e=mgr.Var("e");


    BDD final=(a == b) & (c == d) | e;
    a.clear();
    b.clear();
    c.clear();
    d.clear();
    e.clear();


    mgr.DumpTikZ(std::cout);
}

void test3()
{
    miniBDD_mgr mgr;

    BDD final=mgr.Var("x") & mgr.Var("y");

    //mgr.DumpDot(std::cout);
    //mgr.DumpTikZ(std::cout);
    mgr.DumpTable(std::cout);
}

void test4()
{
    miniBDD_mgr mgr;

    BDD a=mgr.Var("a");
    BDD b=mgr.Var("b");
    BDD c=mgr.Var("c");

    BDD f=!a;

    a.clear();
    b.clear();
    c.clear();


    //mgr.DumpDot(std::cout);
    mgr.DumpTable(std::cout);
    mgr.DumpTikZ(std::cout);
}

void test5()
{
    miniBDD_mgr mgr;

    BDD a=mgr.Var("a");
    BDD b=mgr.Var("b");
    BDD c=mgr.Var("c");

    BDD f=a&!b|!(a|c&b);

    a.clear();
    b.clear();
    c.clear();


    mgr.DumpTable(std::cout);
    mgr.DumpNumber(std::cout);
}

void test6()
{
    miniBDD_mgr mgr;

    BDD a=mgr.Var("a");
    BDD b=mgr.Var("b");
    BDD c=mgr.Var("c");
    BDD d=mgr.Var("d");

    BDD f=(a|!b)&(a&!c|!d)|!a&b&c&d;

    a.clear();
    b.clear();
    c.clear();
    d.clear();


    mgr.DumpTable(std::cout);
    mgr.DumpNumber(std::cout);
}

void test7() //初期サンプル
{
    miniBDD_mgr mgr;

    BDD a=mgr.Var("a");
    BDD b=mgr.Var("b");
    BDD c=mgr.Var("c");

    BDD f=!(a|b)|!(a|c);

    a.clear();
    b.clear();
    c.clear();


    mgr.DumpTable(std::cout);
    mgr.DumpNumber(std::cout);

}

void test8()
{
    miniBDD_mgr mgr;

    BDD a=mgr.Var("a");
    BDD b=mgr.Var("b");
    BDD c=mgr.Var("c");
    BDD d=mgr.Var("d");

    BDD f=a&!b|b&!c&!d|!a&b&c&d;

    a.clear();
    b.clear();
    c.clear();
    d.clear();


    mgr.DumpTable(std::cout);
    mgr.DumpNumber(std::cout);
}

void test9()
{
    miniBDD_mgr mgr;

    BDD a=mgr.Var("a");
    BDD b=mgr.Var("b");
    BDD c=mgr.Var("c");
    BDD d=mgr.Var("d");

    BDD f=a&b&c&d|!a&!b&!c&!d;
    a.clear();
    b.clear();
    c.clear();
    d.clear();

    mgr.DumpTable(std::cout);
    mgr.DumpNumber(std::cout);
}

void test10()
{
    miniBDD_mgr mgr;

    BDD a=mgr.Var("a");
    BDD b=mgr.Var("b");
    BDD c=mgr.Var("c");
    //BDD d=mgr.Var("d");

    BDD f=!a&!b&c|!a&b&!c|a&b&c|a&!b&!c;
    a.clear();
    b.clear();
    c.clear();
    //d.clear();

    mgr.DumpTable(std::cout);
    mgr.DumpNumber(std::cout);
}

void test11()
{
    miniBDD_mgr mgr;

    BDD a=mgr.Var("a");
    BDD b=mgr.Var("b");
    BDD c=mgr.Var("c");
    //BDD d=mgr.Var("d");

    BDD f=a&b|b&c|a&c;
    a.clear();
    b.clear();
    c.clear();
    //d.clear();

    mgr.DumpTable(std::cout);
    mgr.DumpNumber(std::cout);
}

void fourbitx_8()
{
    miniBDD_mgr mgr;
    BDD A=mgr.Var("A");
    BDD B=mgr.Var("B");
    BDD C=mgr.Var("C");
    BDD D=mgr.Var("D");
    BDD E=mgr.Var("E");
    BDD F=mgr.Var("F");
    BDD G=mgr.Var("G");
    BDD H=mgr.Var("H");
    BDD f=!A&D&H | A&!B&C&D&E&!F&!G | B&D&H | C&D&H | D&E&!H | D&F&H | D&G&H;
    A.clear();
    B.clear();
    C.clear();
    D.clear();
    E.clear();
    F.clear();
    G.clear();
    H.clear();
    mgr.DumpTable(std::cout);
    mgr.DumpNumber(std::cout);
}

void fourbitx_7()
{
    miniBDD_mgr mgr;
    BDD A=mgr.Var("A");
    BDD B=mgr.Var("B");
    BDD C=mgr.Var("C");
    BDD D=mgr.Var("D");
    BDD E=mgr.Var("E");
    BDD F=mgr.Var("F");
    BDD G=mgr.Var("G");
    BDD H=mgr.Var("H");
    BDD f=!A&D&G&!H | A&!B&!C&D&E&!F&H | B&D&G&!H | !C&D&G | C&!D&H | C&!G&H | D&E&G&!H | D&F&G&!H;
    A.clear();
    B.clear();
    C.clear();
    D.clear();
    E.clear();
    F.clear();
    G.clear();
    H.clear();
    mgr.DumpTable(std::cout);
    mgr.DumpNumber(std::cout);
}

void fourbitx_6()
{
    miniBDD_mgr mgr;
    BDD A=mgr.Var("A");
    BDD B=mgr.Var("B");
    BDD C=mgr.Var("C");
    BDD D=mgr.Var("D");
    BDD E=mgr.Var("E");
    BDD F=mgr.Var("F");
    BDD G=mgr.Var("G");
    BDD H=mgr.Var("H");
    BDD f=!B&C&D&F | !B&C&!D&G | !B&C&F&G&H | !B&D&F&!G | B&!C&!D&H | B&!C&!F&H | B&C&D&!F&G | B&!D&!G&H | B&!F&!G&H | !C&D&F&!H | C&!D&G&!H | C&!F&G&!H | D&F&!G&!H;
    A.clear();
    B.clear();
    C.clear();
    D.clear();
    E.clear();
    F.clear();
    G.clear();
    H.clear();
    mgr.DumpTable(std::cout);
    mgr.DumpNumber(std::cout);
}

void fourbitx_5()
{
    miniBDD_mgr mgr;
    BDD A=mgr.Var("A");
    BDD B=mgr.Var("B");
    BDD C=mgr.Var("C");
    BDD D=mgr.Var("D");
    BDD E=mgr.Var("E");
    BDD F=mgr.Var("F");
    BDD G=mgr.Var("G");
    BDD H=mgr.Var("H");
    BDD f=!A&!B&!C&D&E | !A&!B&C&!D&F | !A&!B&C&D&!E&!F&G&H | !A&!B&C&!E&F&!G | !A&!B&D&E&F&G | !A&B&!C&!D&G | !A&B&!C&D&!E&F&!G&H | !A&B&!C&!E&!F&G | !A&B&C&D&E&H | !A&B&D&F&G&H | !A&B&E&F&G&H | !A&C&!D&F&!G | !A&D&E&!F&!G | A&!B&!C&!D&H | A&!B&!C&!E&H | A&!B&C&D&E&F&!G&H | A&!B&C&E&!F&G&H | A&!B&!D&!F&H | A&B&!C&D&E&!F&G&H | A&B&!C&E&F&!G&H | A&B&C&!D&!F&G | A&B&C&D&!E&F | A&B&C&!E&!F&H | A&!C&!D&!G&H | A&!D&!F&!G&H | A&D&!E&F&G&H | A&!E&!F&!G&H | !B&!C&D&E&!H | !B&C&!D&F&!H | !B&C&E&F&G&!H | !B&D&E&!F&!H | B&!C&!D&G&!H | B&!C&!E&G&!H | B&C&D&!E&F&!H | B&!D&!F&G&!H | B&!E&!F&G&!H | !C&D&E&!G&!H | C&!D&F&!G&!H | C&!E&F&!G&!H | D&E&!F&!G&!H;
    A.clear();
    B.clear();
    C.clear();
    D.clear();
    E.clear();
    F.clear();
    G.clear();
    H.clear();
    mgr.DumpTable(std::cout);
    mgr.DumpNumber(std::cout);
}

void fourbitx_4()
{
    miniBDD_mgr mgr;
    BDD A=mgr.Var("A");
    BDD B=mgr.Var("B");
    BDD C=mgr.Var("C");
    BDD D=mgr.Var("D");
    BDD E=mgr.Var("E");
    BDD F=mgr.Var("F");
    BDD G=mgr.Var("G");
    BDD H=mgr.Var("H");
    BDD f=A&!B&C&!D&E | A&!B&C&E&!F&!G | A&B&!C&!D&F | A&B&!C&D&E&!F&G | A&B&!C&!E&F&!G | A&B&!C&!E&F&!H | A&B&C&D&!E&G | A&B&!D&!E&F&!G | A&B&D&E&F&!G&!H | A&C&!D&E&!F&!G | A&C&!D&E&F&G | A&C&D&!E&F&G | A&C&E&!F&!G&!H | A&!B&!C&!D&G&!H | A&!B&!C&D&E&!G | A&!B&!C&!E&G&!H | A&!B&!C&!G&H | A&!B&C&!E&F&!G&!H | A&!B&!D&!E&!F&G&!H | A&B&!C&!D&!F&H | A&B&C&D&!F&H | A&B&C&!E&H | A&B&!D&E&F&!H | A&B&!E&!F&H | A&!C&D&E&!F&!G&!H | A&!C&!E&F&G&!H | A&C&!D&!E&F&H | A&D&E&F&!G&H | A&!E&!F&!G&H | B&C&!D&E&!F&G&H | B&C&D&E&!F&G&!H | B&!C&!D&F&!H | B&C&!D&E&!F&!H | B&C&!E&!F&G&H | B&D&E&!F&G&H | C&!D&E&!F&!G&!H | C&!D&E&F&G&!H | C&D&!E&F&G&H;
    A.clear();
    B.clear();
    C.clear();
    D.clear();
    E.clear();
    F.clear();
    G.clear();
    H.clear();
    mgr.DumpTable(std::cout);
    mgr.DumpNumber(std::cout);
}

void fourbitx_3()
{
    miniBDD_mgr mgr;
    BDD A=mgr.Var("A");
    BDD B=mgr.Var("B");
    BDD C=mgr.Var("C");
    BDD D=mgr.Var("D");
    BDD E=mgr.Var("E");
    BDD F=mgr.Var("F");
    BDD G=mgr.Var("G");
    BDD H=mgr.Var("H");
    BDD f=A&!B&C&D&E&F | A&B&!C&!D&E | A&B&!C&E&!F | A&B&!C&E&!G&!H | A&B&C&D&!E&F&H | A&B&C&!E&F&G | A&B&!D&E&!F&!H | A&B&D&!E&F&G&H | A&B&E&!F&!G | A&C&D&E&F&G | A&!B&!C&!D&F | A&!B&!C&!E&F | A&!B&!C&F&!H | A&!B&C&D&E&!F&H | A&!B&C&E&!F&G | A&!B&!D&!E&F&!H | A&!B&D&E&!F&G&H | A&!B&!E&F&!G | A&B&!C&D&E&!G&H | A&B&C&!D&E&!G | A&B&C&F&!G&!H | A&B&!E&!F&G&H | A&!C&!D&!E&F&!G | A&!C&D&E&F&!G&H | A&!C&E&F&G&!H | A&C&!D&F&!G&!H | A&C&D&!F&G&H | A&!E&F&!G&!H | B&C&D&E&G&H | B&!C&!D&E&F&G | B&!C&!D&E&G&!H | B&C&!E&F&G&H | B&!D&E&!F&!G | B&E&!F&!G&!H | C&D&E&F&G&H;
    A.clear();
    B.clear();
    C.clear();
    D.clear();
    E.clear();
    F.clear();
    G.clear();
    H.clear();
    mgr.DumpTable(std::cout);
    mgr.DumpNumber(std::cout);
}

void fourbitx_2()
{
    miniBDD_mgr mgr;
    BDD A=mgr.Var("A");
    BDD B=mgr.Var("B");
    BDD C=mgr.Var("C");
    BDD D=mgr.Var("D");
    BDD E=mgr.Var("E");
    BDD F=mgr.Var("F");
    BDD G=mgr.Var("G");
    BDD H=mgr.Var("H");
    BDD f=A&B&C&D&E&G | A&B&C&E&F | A&B&C&E&G&H | A&B&D&E&F&G | A&B&D&E&F&H | A&!B&!C&!D&E | A&!B&!C&E&!G | A&!B&!C&E&!H | A&!B&!D&E&!G&!H | A&!B&E&!F | A&B&C&!E&F&H | A&B&D&!E&F&H | A&B&!E&F&G | A&!C&!D&E&!F&!H | A&!C&E&!F&!G | A&C&D&!E&F&G | A&C&!E&F&G&H | A&!D&E&!F&!G | A&E&!F&!G&!H | B&C&D&E&F&H | B&C&E&F&G | B&D&E&F&G&H;
    A.clear();
    B.clear();
    C.clear();
    D.clear();
    E.clear();
    F.clear();
    G.clear();
    H.clear();
    mgr.DumpTable(std::cout);
    mgr.DumpNumber(std::cout);
}

void fourbitx_1()
{
    miniBDD_mgr mgr;
    BDD A=mgr.Var("A");
    BDD B=mgr.Var("B");
    BDD C=mgr.Var("C");
    BDD D=mgr.Var("D");
    BDD E=mgr.Var("E");
    BDD F=mgr.Var("F");
    BDD G=mgr.Var("G");
    BDD H=mgr.Var("H");
    BDD f=A&B&C&D&E&H | A&B&C&E&G | A&B&D&E&G | A&B&E&F | A&B&E&G&H | A&C&D&E&F | A&C&E&F&G | A&C&E&F&H | A&D&E&F&G&H;
    A.clear();
    B.clear();
    C.clear();
    D.clear();
    E.clear();
    F.clear();
    G.clear();
    H.clear();
    mgr.DumpTable(std::cout);
    mgr.DumpNumber(std::cout);
}

void threebitx_6()
{
    miniBDD_mgr mgr;
    BDD A=mgr.Var("A");
    BDD B=mgr.Var("B");
    BDD C=mgr.Var("C");
    BDD D=mgr.Var("D");
    BDD E=mgr.Var("E");
    BDD F=mgr.Var("F");
    BDD f=!A&B&!C&D | !A&B&C&!D&E&F | !A&B&D&!E | A&!B&!C&E | A&!B&C&D&E&F | A&!B&!D&E | A&!B&E&!F | A&B&C&D&!F | A&!D&E&!F | B&!C&D&!E | B&!C&D&F ;
    A.clear();
    B.clear();
    C.clear();
    D.clear();
    E.clear();
    F.clear();

    mgr.DumpTable(std::cout);
    mgr.DumpNumber(std::cout);
}

int main()
{
    fourbitx_4();
}