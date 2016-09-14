#include "d:\work\bc\include\mcadincl.h"
#include <math.h>
#include <winuser.h>
#include <stdio.h>
#include <alloc.h>
#include <string.h>
#include <io.h>






    LRESULT FDNAFunction(   COMPLEXARRAY * const E_D,
    								 COMPLEXSCALAR * const Pd_Pm,
                            COMPLEXSCALAR * const A,							// Определяет Коэффициент усиления и КНД антенны
                            COMPLEXARRAY * const delta,						// E_D[0] - если Pd_Pm==0 - коэффициент усиления (E) передающей антенны,
                            COMPLEXSCALAR * const l,							// 			если Pd_Pm==1 - КНД приемной антенны
                            COMPLEXARRAY * const lambda,               // E_D[1] - тип антенны буквы закодированы числами
                            COMPLEXSCALAR * const lambda_0,             // Pd_Pm - 0 - производить расчет коэффициента усиления (E) передающей антенны,
                            COMPLEXSCALAR * const H);                   //         1 - производить расчет КНД приемной антенны
                                                                        // A - Азимут на корреспондента относительно главного лепестка,
																								// delta - Массив углов возвышения излучения (прихода) волны в вертикальной
                                                                        // плоскости из 18 элементов,
                                                                        // l - Длина плеча антенны,
                                                                        // lambda - Массив рабочих длин волн из 18 элементов,
                                                                        // lambda_0 - Резонансная длина волны (для настроенных антенн)
                                                                        // H - Высота установки антенны


    double linterp(int N,                   // Линейная интерполяция функции по точкам. За внешними пределами вектора X выдает значения крайних точек
    						double* X,               // X - вектор абсцисс,
                     double* Y,               // Y - вектор ординат,
                     double x);               // x - значение абсциссы точки, в которой производится вычисление
                                              // N - количество точек

    char *FReadDNA_KVAFunction(double *DNA_KVA,
    										int Pd_Pm);			 // Читает файл диаграммы направленности коротковолновой антенны и построчно копирует его в массив
                                 						 //DNA_KVA.
                                                    // Pd_Pm - 0 - передающая антенна, 1 - приемная антенна.
                                                    // Возвращает указатель на строку, содержащую тип антенны
                                                    //В случае ошибки возвращает NULL

FUNCTIONINFO    FDNA =
    { 
    "FDNA",                          // Name by which mathcad will recognize the function
    "Pd_Pm,A,delta,l,lambda,lambda_0,H",    // FDNA will be called as FDNA(DNA,delta,l,lambda,lambda_0,H)
    "Определяет Коэффициент усиления и КНД антенны",      // description of FDNA(DNA,A,delta,l,lambda,lambda_0,H)
    (LPCFUNCTION)FDNAFunction,       // pointer to the executible code
    COMPLEX_ARRAY,                     // the return type is also a complex array
    7,                                  //  the function takes on 7 arguments
    { COMPLEX_SCALAR,COMPLEX_SCALAR, COMPLEX_ARRAY,COMPLEX_SCALAR,COMPLEX_ARRAY,COMPLEX_SCALAR,COMPLEX_SCALAR}
    };

    double linterp(int N,                   // Линейная интерполяция функции по точкам.  За внешними пределами вектора X выдает значения крайних точек
    						double* X,               // X - вектор абсцисс,
                     double* Y,               // Y - вектор ординат,
                     double x)               // x - значение абсциссы точки, в которой производится вычисление
                                              // N - количество точек
    {
    	double y;
      int i, x_vozr;

      if(X[0]<X[1])
      	x_vozr=1;//аргумент возрастает
      else
      	x_vozr=0;//аргумент убывает

      if(x_vozr==1)//возрастающий порядок аргумента
      {
         if(x<=X[0])//если аргумент искомого значения левее 0-й точки
         {
         	y=Y[0];
         }
         else if(x>=X[N-1])//иначе если аргумент искомого значения правее последней точки
         {
         	y=Y[N-1];
         }
         else//Иначе, если аргумент искомого значения правее 0-й и левее последней точки
      		for(i=0;i<=N-2;i++)
      		{
         		if(X[i]<x && x<=X[i+1])
               {
                  y=Y[i]+(Y[i+1]-Y[i])/(X[i+1]-X[i])*(x-X[i]);
         			break;
               }
      		}
      }
      else//иначе, если убывающий порядок аргумента
      {
         if(x<=X[N-1])//если аргумент искомого значения левее (меньше аргумента) последней точки
         {
         	y=Y[N-1];
         }
         else if(x>=X[0])//иначе если аргумент искомого значения правее (наибольшее значение) 0-й точки
         {
         	y=Y[0];
         }
         else//Иначе, если аргумент искомого значения левее 0-й и правее последней точки
      		for(i=0;i<=N-2;i++)
      		{
         		if(X[i]>x && x>=X[i+1])
               {
                  y=Y[i+1]+(Y[i]-Y[i+1])/(X[i]-X[i+1])*(x-X[i+1]);
         			break;
               }
      		}
      }

      return y;
    }












    
    LRESULT FDNAFunction(   COMPLEXARRAY * const E_D,
    								 COMPLEXSCALAR * const Pd_Pm,
                            COMPLEXSCALAR * const A,							// Определяет Коэффициент усиления и КНД антенны
                            COMPLEXARRAY * const delta,						// E_D[0] - если Pd_Pm==0 - коэффициент усиления (E) передающей антенны,
                            COMPLEXSCALAR * const l,							// 			если Pd_Pm==1 - КНД приемной антенны
                            COMPLEXARRAY * const lambda,               // E_D[1] - тип антенны буквы закодированы числами
                            COMPLEXSCALAR * const lambda_0,             // Pd_Pm - 0 - производить расчет коэффициента усиления (E) передающей антенны,
                            COMPLEXSCALAR * const H)                   //         1 - производить расчет КНД приемной антенны
                                                                        // A - Азимут на корреспондента относительно главного лепестка,
																								// delta - Массив углов возвышения излучения (прихода) волны в вертикальной
                                                                        // плоскости из 18 элементов,
                                                                        // l - Длина плеча антенны,
                                                                        // lambda - Массив рабочих длин волн из 18 элементов,
                                                                        // lambda_0 - Резонансная длина волны (для настроенных антенн)
                                                                        // H - Высота установки антенны
    {
        int i,j,k,o,*in,i_e,i_d,*i_e_em,n_e_em=0,n_pod,*n_lines_pod,/* *nomer_first_lin_sem*/*i_lin,*inT,nt,**index_lin_podsem,
        		Index_N_semeystv=2,*Index_N_parametrov,*Index_N_liniy,**Index_N_tochek, Pd_Pm_int, pod_exist;
        double ****XY_DNA,***P_DNA,phi[18]/*заложено, но не реализовано*/,Abscissa[18][6]/*18-количество расчетных частот, 6-количество типов абсцисс*/,
        			Parameter[18][8]/*18-количество расчетных частот, 8-количество типов параметров*/,
               **par_sem,*par_podsem,*e,*e1,E_m[18],D_m[18],E_Em[18]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},*e_em[18];
        double DNA[4096]; //Массив данных диаграмм направленности (описание в файле "Формат файла KVA.txt"),
        char *DNA_KVA_name;//Указатель на тип антенны

        Pd_Pm_int=(int)Pd_Pm->real;

       DNA_KVA_name=FReadDNA_KVAFunction(DNA,
    										Pd_Pm_int);			 // Читает файл диаграммы направленности коротковолновой антенны и построчно копирует его в массив
                               						 //DNA.
                                                    // Pd_Pm_int - 0 - передающая антенна, 1 - приемная антенна.
                                                    // Возвращает указатель на строку, содержащую тип антенны
                                                    //В случае ошибки возвращает NULL

         if(DNA_KVA_name==NULL)
         		return 4;  //Ошибка чтения файла диаграммы направленности антенны

        // allocate space for the return array E_D
        if ( !MathcadArrayAllocate( E_D,  // allocate space for E_D
                    (strlen(DNA_KVA_name)+1<18)?18:strlen(DNA_KVA_name)+1,    //with 18 or strlen(DNA_KVA_name)+1 rows
                    2,    //and 2 cols
                    TRUE,   //  allocate the real part
                    FALSE    //  don't allocate the imaginary part
                    ) )
            return 2;           // if allocation is insufficient
                                // return the error code

        in=(int*)calloc(DNA[Index_N_semeystv], sizeof(int));//массив элементов int по количеству семейств (индексы количества типов параметров каждого семейства)

        in[0]=5;//5- индекс количества типов параметров 0-го семейства
        for(i=1;i<=DNA[Index_N_semeystv]-1;i++)//перебрать семейства кривых начиная с 1-го (i - номер семейства)
        {
         inT=(int*)calloc(DNA[in[i-1]+DNA[in[i-1]]+1], sizeof(int));//массив элементов int по количеству линий в i-1-м семействе
        	inT[0]=in[i-1]+2*DNA[in[i-1]]+2;//индекс количества точек в 0-й линии в i-1-м семействе
         nt=DNA[inT[0]];//количество точек в 0-й линии в i-1-м семействе
         if(DNA[in[i-1]+DNA[in[i-1]]+1]>1)//если количество линий в i-1-м семействе больше 1
         	for(j=1;j<=DNA[in[i-1]+DNA[in[i-1]]+1]-1;j++)//перебрать линии i-1-го семейства начиная с 1-й
         	{
         		inT[j]=inT[j-1]+2*DNA[inT[j-1]]+DNA[in[i-1]]+1;//индекс количества точек в остальных линиях в i-1-м семействе
            	nt=nt+DNA[inT[j]];//суммарное количество точек в i-1-м семействе
         	}
         free(inT);
         //Индекс количества типов параметров i-го семейства
         in[i]=in[i-1]+DNA[in[i-1]]*(DNA[in[i-1]+DNA[in[i-1]]+1]+1)+4+DNA[in[i-1]+DNA[in[i-1]]+1]+2*nt;
        }
        Index_N_parametrov=in;//массив индексов количества типов параметров
        Index_N_liniy=(int*)calloc(DNA[Index_N_semeystv], sizeof(int));//массив индексов количества линий в семействах
        for(i=0;i<=DNA[Index_N_semeystv]-1;i++)//перебрать все семейства кривых (i - номер семейства)
        	Index_N_liniy[i]=Index_N_parametrov[i]+DNA[Index_N_parametrov[i]]+1;//индекс количества линий в i-том семействе
        Index_N_tochek=(int**)calloc(DNA[Index_N_semeystv], sizeof(int*));//массив указателей на int по количеству семейств

        for(i=0;i<=DNA[Index_N_semeystv]-1;i++)//перебрать все семейства кривых (i - номер семейства)
        {
        	Index_N_tochek[i]=(int*)calloc(DNA[Index_N_liniy[i]], sizeof(int));//массив указателей по количеству линий в i-том семействе
         Index_N_tochek[i][0]=Index_N_liniy[i]+DNA[Index_N_parametrov[i]]+1;//индекс количества точек в 0-й линии i-го семейства
         if(DNA[Index_N_liniy[i]]>1)//если количество линий в i-том семействе больше 1-й
         	for(j=1;j<=DNA[Index_N_liniy[i]]-1;j++)//перебрать линии i-го семейства начиная с 1-й (j - номер линии)
         		Index_N_tochek[i][j]=Index_N_tochek[i][j-1]+2*DNA[Index_N_tochek[i][j-1]]+DNA[Index_N_parametrov[i]]+1;//индекс количества точек в j-й линии i-го семейства
        }
        XY_DNA=(double****)calloc(DNA[Index_N_semeystv], sizeof(double***));//массив указателей на double*** по количеству семейств
        for(i=0;i<=DNA[Index_N_semeystv]-1;i++)//перебрать все семейства (i - номер семейства)
        {
        	XY_DNA[i]=(double***)calloc(DNA[Index_N_liniy[i]], sizeof(double**));//массив указателей на double** по количеству линий в i-том семействе
         for(j=0;j<=DNA[Index_N_liniy[i]]-1;j++)//перебрать все линии i-го семейства (j - номер линии)
         {
         	XY_DNA[i][j]=(double**)calloc(2, sizeof(double*));//массив указателей на double* из 2-х элементов (указатель на абсциссы и указатель на ординаты j-той линии i-го семейства)
            XY_DNA[i][j][0]=(double*)calloc(DNA[Index_N_tochek[i][j]], sizeof(double));//массив значений абсцисс по количеству точек в j-той линии i-го семейства
            XY_DNA[i][j][1]=(double*)calloc(DNA[Index_N_tochek[i][j]], sizeof(double));//массив значений ординат по количеству точек в j-той линии i-го семейства
            for(k=0;k<=DNA[Index_N_tochek[i][j]]-1;k++)//перебрать все точки j-той линии i-го семейства (k - номер точки)
            {
            	XY_DNA[i][j][0][k]=DNA[Index_N_tochek[i][j]+1+2*k];//значение абсциссы k-той точки, j-той линии, i-того семейства
               XY_DNA[i][j][1][k]=DNA[Index_N_tochek[i][j]+2+2*k];//значение ординаты k-той точки, j-той линии, i-того семейства
            }
         }
        }
        P_DNA=(double***)calloc(DNA[Index_N_semeystv], sizeof(double**));//массив указателей на double** по количеству семейств
        for(i=0;i<=DNA[Index_N_semeystv]-1;i++)//перебрать все семейства (i - номер семейства)
        {
        	P_DNA[i]=(double**)calloc(DNA[Index_N_liniy[i]], sizeof(double*));//массив указателей на double* по количеству линий в i-том семействе
         for(j=0;j<=DNA[Index_N_liniy[i]]-1;j++)//перебрать все линии i-го семейства
         {
         	P_DNA[i][j]=(double*)calloc(DNA[Index_N_parametrov[i]], sizeof(double));//массив значений параметров для j-той линии в i-том семействе
            if(DNA[Index_N_parametrov[i]]>0)//если в i-том семействе больше одной линии
            	for(k=0;k<=DNA[Index_N_parametrov[i]]-1;k++)//перебрать все параметры i-го семейства (k - номер типа параметра)
            		P_DNA[i][j][k]=DNA[Index_N_tochek[i][j]-DNA[Index_N_parametrov[i]]+k];//значение k-того параметра, j-той линии, i-того семейства
         }
        }

        i_e_em=(int*)calloc(DNA[Index_N_semeystv],sizeof(int));//массив номеров семейств типа E/Емакс
        for(i=0;i<=DNA[Index_N_semeystv]-1;i++)//перебрать все семейства (i - номер семейства)
        {
         switch ((short)DNA[Index_N_parametrov[i]-1])//определение типа оси ординат для i-го семейства (по сути - типа семейства)
         {
         	case 0: i_d=i; break;//если тип оси ординат для i-го семейства - D(КНД), то i - номер семейства типа D
         	case 1: i_e=i; break;//если тип оси ординат для i-го семейства - Е(Коэффициент усиления), то i - номер семейства типа E
         	case 2: i_e_em[n_e_em]=i; n_e_em++; break;//если тип оси ординат для i-го семейства - Е/Емакс(нормированная ДНА в горизонтальной плоскостия),
            														//то i - номер n_e_em-того семейства типа E/Емакс
         	case 3:  for(k=0; k<18; k++)//для всех расчетных частот
            				phi[k]=linterp(DNA[Index_N_tochek[i][0]],XY_DNA[i][0][0],XY_DNA[i][0][1],300./lambda->hReal[0][k]); break;//если тип оси ординат
            														//для i-го семейства - Phi (град., угол фазирования антенн типа БС),
            														//то вычисляем значение Phi
         	default:;
         }
        }
        for(i=0; i<18; i++)//для всех расчетных частот
        {
        	Abscissa[i][0]=lambda->hReal[0][i]/lambda_0->real;//значение отношения расчетной длины волны к оптимальной (резонансной) длине волны антенны
        	Abscissa[i][1]=lambda->hReal[0][i];//м, значение расчетной длины волны
        	Abscissa[i][2]=A->real;//град., угол направления в горизонтальной плоскости отн-но главного лепестка	ДН
        	Abscissa[i][3]=delta->hReal[0][i];//град., угол прихода (излучения) волны в вертикальной плоскости отн-но касательной к поверхности земли
        	Abscissa[i][4]=l->real/lambda->hReal[0][i];//отношение длины плеча антенны к длине волны
        	Abscissa[i][5]=300/lambda->hReal[0][i];//МГц, расчетная частота
        	Parameter[i][0]=delta->hReal[0][i];//град., угол прихода (излучения) волны в вертикальной плоскости отн-но земли
        	Parameter[i][1]=l->real;//м, длина плеча антенны
        	Parameter[i][2]=lambda->hReal[0][i]/lambda_0->real;//значение отношения расчетной длины волны к оптимальной (резонансной) длине волны антенны
        	Parameter[i][3]=lambda->hReal[0][i];//м, расчетная длина волны
        	Parameter[i][4]=phi[i];//град., угол фазирования антенн типа БС
        	Parameter[i][5]=H->real/lambda_0->real;//отношение высоты установки антенны к резонансной длине волны
        	Parameter[i][6]=l->real/lambda->hReal[0][i];//отношение длины плеча антенны к длине волны
        	Parameter[i][7]=H->real/lambda->hReal[0][i];//отношение высоты установки антенны к длине волны
        }

        //Определение значения коэффициента усиления передающей антенны (Е)
        if(Pd_Pm_int==0)//если передающая антенна
        {
        	if(DNA[Index_N_parametrov[i_e]]>1)//если количество типов параметров в семействе Е больше 1, т.е. 2
        	{
        		n_pod=1;//количество подсемейств (подсемейства группируются по равенству второго (с индексом 1) типа параметра семейства).
            		  //начальное значение - 1 подсемейство
	        	n_lines_pod=(int*)calloc(DNA[Index_N_liniy[i_e]],sizeof(int));//массив значений количества линий в подсемействе семейства Е
   	     	n_lines_pod[0]=1;//количество линий в 0-м подсемействе семейства Е (начальное значение - 1 линия в 0-м подсемействе)

            par_podsem=(double*)calloc(DNA[Index_N_liniy[i_e]],sizeof(double));//массив значений параметров второго типа (с индексом 1)
            																						 //подсемейства семейства Е
        		par_podsem[0]=P_DNA[i_e][0][1];//значение параметра второго типа (с индексом 1) 0-го подсемейства семейства Е
            										 //(равен параметру второго типа 0-й линии семейства E)
	         /*Код, замененный на новый
      	  	nomer_first_lin_sem=(int*)calloc(DNA[Index_N_liniy[i_e]],sizeof(int));//указатель на номер нулевой линии в семействе E
	         nomer_first_lin_sem[0]=0;//номер нулевой линии в семействе
	         */

            index_lin_podsem=(int**)calloc(DNA[Index_N_liniy[i_e]],sizeof(int*));//массив индексов линий подсемейств семейства Е,
            																									 //соответствующих индексам в массивах P_DNA и XY_DNA
            for(i=0; i<=DNA[Index_N_liniy[i_e]]-1;i++)//выделить память по количеству линий семейства Е
            	index_lin_podsem[i]=(int*)calloc(DNA[Index_N_liniy[i_e]],sizeof(int));//массив индексов линий i-го подсемейства семейства Е,
            																									 //соответствующих индексам линиий в массивах P_DNA и XY_DNA
            index_lin_podsem[0][0]=0;//индекс 0-ой линий (второй индекс) 0-го подсемейства (первый индекс) семейства Е,
            								 //соответствующий индексам линий в массивах P_DNA и XY_DNA

            for(i=1;i<=DNA[Index_N_liniy[i_e]]-1;i++)//перебрать линии семейства Е начиная с 1-й (информация о 0-й линии уже занесена в массивы)
	         {
               /*Код, замененный на новый
	            if(P_DNA[i_e][i][1]!=P_DNA[i_e][i-1][1])//если значение 1-го параметра i-й линии семейства Е не равно
               													 //1-му параметру i-1-й линии семейства Е
	         	{
	            	par_podsem[n_pod]=P_DNA[i_e][i][1];//значение параметра n_pod-го подсемейства семейства Е
	               nomer_first_lin_sem[n_pod]=i;//номер первой линии n_pod-го подсемейства семейства Е
	               n_lines_pod[n_pod-1]=1;//количество линий n_pod-1-го параметра подсемейства семейства Е
	               n_pod++;
	            }
	            else
	            	n_lines_pod[n_pod-1]++;
               */

               pod_exist=0;//признак, что подсемейство со значением 1-го параметра текущей (i-той) линии семейства Е существует
               				//0 - не существует, 1 - существует
               for(j=0; j<n_pod; j++)//Перебрать все найденные подсемейства семейства Е
               {
               	if(P_DNA[i_e][i][1]==P_DNA[i_e][index_lin_podsem[j][0]][1])//Если значение 1-го параметра i-ой линии семейства Е == значению
                  																						//1-го параметра 0-ой линии j-го подсемейства (т.е. подсемейство
                                                                                    //с таким значением 1-го параметра уже есть)
                  {
                     index_lin_podsem[j][n_lines_pod[j]]=i;//индекс текущей линии j-го подсемейства,соответствующий индексам в массивах P_DNA и XY_DNA
                     n_lines_pod[j]++;//увеличение количества линий в j-том подсемействе
                     pod_exist=1;//подсемейство со значением 1-го параметра текущей (i-той) линии семейства Е существует
                  }
               }
               if(pod_exist==0)//если подсемейства с таким значением первого параметра еще нет
               {
                  par_podsem[n_pod]=P_DNA[i_e][i][1];//значение параметра второго типа (с индексом 1) n_pod-го подсемейства семейства Е
            										 //(равен параметру второго типа i-й линии семейства E)
               	index_lin_podsem[n_pod][0]=i;//i - индекс 0-й линии нового (n_pod-го) подсемейства,соответствующий индексам в массивах P_DNA и XY_DNA
                  n_lines_pod[n_pod]=1;//количества линий в новом подсемействе - 1
                  n_pod++;//увеличение количества подсемейств
               }
	         }
 	      	/*k=0;*/
	         par_sem=(double**)calloc(n_pod,sizeof(double*));//указатель на массив значений параметров 0-го типа внутри подсемейства
	         for(i=0;i<=n_pod-1;i++)//перебрать все подсемейства
	         {
	            par_sem[i]=(double*)calloc(n_lines_pod[i],sizeof(double));//массив значений параметров 0-го типа внутри подсемейства
	         	for(j=0;j<=n_lines_pod[i]-1;j++)//перебрать все линии i-го подсемейства
 	           	{
	            	par_sem[i][j]=P_DNA[i_e][index_lin_podsem[i][j]][0];
              	}
	         }
            for(k=0; k<18; k++)//по всем рабочим длинам волн
            {
            	e=(double*)calloc(n_pod,sizeof(double));//указатель на массив значений Е для каждого подсемейства
	         	for(i=0;i<=n_pod-1;i++)//по всем подсемействам
	         	{
	         		e1=(double*)calloc(n_lines_pod[i],sizeof(double));//указатель на массив значений Е для каждой линии подсемейства
	            	for(j=0;j<=n_lines_pod[i]-1;j++)//по всем линиям подсемейства
                     //Интерполяция значений Е для каждой линии внутри i-го подсемейства по значениям абсцисс для каждой рабочей длины волны
	            		e1[j]=linterp(DNA[Index_N_tochek[i_e][index_lin_podsem[i][j]]],XY_DNA[i_e][index_lin_podsem[i][j]][0],XY_DNA[i_e][index_lin_podsem[i][j]][1],Abscissa[k][DNA[Index_N_parametrov[i_e]-2]]);
                  //Интерполяция значений Е внутри i-го подсемейства по первому (с индексом 0) параметру
	            	e[i]=linterp(n_lines_pod[i],par_sem[i],e1,Parameter[k][DNA[Index_N_parametrov[i_e]+1]]);
	            	free(e1);
	         	}
               //Интерполяция значений Е между подсемействами по второму (с индексом 1) параметру
	         	E_m[k]=linterp(n_pod,par_podsem,e,Parameter[k][DNA[Index_N_parametrov[i_e]+2]]);
	         	free(e);
            }
	         for(i=0;i<=n_pod-1;i++)
	         	free(par_sem[i]);
	         free(par_sem);
	         /*free(nomer_first_lin_sem);*/
	         free(par_podsem);
	         free(n_lines_pod);
            for(i=0; i<=DNA[Index_N_liniy[i_e]]-1;i++)//освободить память по количеству линий семейства Е начиная
            	free(index_lin_podsem[i]);
            free(index_lin_podsem);
         }
	      else if(DNA[Index_N_parametrov[i_e]]==1)//если количество параметров в семействе Е - 1
	      {
	         par_sem=(double**)calloc(1,sizeof(double*));
	        	par_sem[0]=(double*)calloc(DNA[Index_N_liniy[i_e]],sizeof(double));
	         for(i=0;i<=DNA[Index_N_liniy[i_e]]-1;i++)
	         	par_sem[0][i]=P_DNA[i_e][i][0];//массив значений 0-го параметря для i-тых линий семейства Е
            for(k=0; k<18; k++)//по всем рабочим длинам волн
            {
	         	e=(double*)calloc(DNA[Index_N_liniy[i_e]],sizeof(double));
	         	for(i=0;i<=DNA[Index_N_liniy[i_e]]-1;i++)
	         		e[i]=linterp(DNA[Index_N_tochek[i_e][i]],XY_DNA[i_e][i][0],XY_DNA[i_e][i][1],Abscissa[k][DNA[Index_N_parametrov[i_e]-2]]);
	         	E_m[k]=linterp(DNA[Index_N_liniy[i_e]],par_sem[0],e,Parameter[k][DNA[Index_N_parametrov[i_e]+1]]);
	         	free(e);
            }
	         free(par_sem[0]);
	         free(par_sem);
	      }
	      else//если количество параметров в семействе Е - 0 (1 линия)
            for(k=0; k<18; k++)//по всем рабочим длинам волн
	      		E_m[k]=linterp(DNA[Index_N_tochek[i_e][0]],XY_DNA[i_e][0][0],XY_DNA[i_e][0][1],Abscissa[k][DNA[Index_N_parametrov[i_e]-2]]);
        }

        //Определение значения КНД приемной антенны (D)
        else//иначе (если приемная антенна)
        {
         if(DNA[Index_N_parametrov[i_d]]>1)//если количество типов параметров в семействе D больше 1, т.е. 2
        	{
        		n_pod=1;//количество подсемейств (подсемейства группируются по равенству второго (с индексом 1) типа параметра семейства).
            		  //начальное значение - 1 подсемейство
	        	n_lines_pod=(int*)calloc(DNA[Index_N_liniy[i_d]],sizeof(int));//массив значений количества линий в подсемействе семейства D
   	     	n_lines_pod[0]=1;//количество линий в 0-м подсемействе семейства D (начальное значение - 1 линия в 0-м подсемействе)

            par_podsem=(double*)calloc(DNA[Index_N_liniy[i_d]],sizeof(double));//массив значений параметров второго типа (с индексом 1)
            																						 //подсемейства семейства D
        		par_podsem[0]=P_DNA[i_d][0][1];//значение параметра второго типа (с индексом 1) 0-го подсемейства семейства D
            										 //(равен параметру второго типа 0-й линии семейства D)
            index_lin_podsem=(int**)calloc(DNA[Index_N_liniy[i_d]],sizeof(int*));//массив индексов линий подсемейств семейства D,
            																									 //соответствующих индексам в массивах P_DNA и XY_DNA
            for(i=0; i<=DNA[Index_N_liniy[i_d]]-1;i++)//выделить память по количеству линий семейства D начиная
            	index_lin_podsem[i]=(int*)calloc(DNA[Index_N_liniy[i_d]],sizeof(int));//массив индексов линий i-го подсемейства семейства D,
            																									 //соответствующих индексам в массивах P_DNA и XY_DNA
            index_lin_podsem[0][0]=0;//индекс 0-ой линий (второй индекс) 0-го подсемейства (первый индекс) семейства D,
            								 //соответствующий индексам линий в массивах P_DNA и XY_DNA
            for(i=1;i<=DNA[Index_N_liniy[i_d]]-1;i++)//перебрать линии семейства D начиная с 1-й (информация о 0-й линии уже занесена в массивы)
	         {
               pod_exist=0;//признак, что подсемейство со значением 1-го параметра текущей (i-той) линии семейства D существует
               				//0 - не существует, 1 - существует
               for(j=0; j<n_pod; j++)//Перебрать все найденные подсемейства семейства D
               {
               	if(P_DNA[i_d][i][1]==P_DNA[i_d][index_lin_podsem[j][0]][1])//Если значение 1-го параметра i-ой линии семейства D == значению
                  																						//1-го параметра 0-ой линии j-го подсемейства (т.е. подсемейство
                                                                                    //с таким значением 1-го параметра уже есть)
                  {
                     index_lin_podsem[j][n_lines_pod[j]]=i;//индекс текущей линии j-го подсемейства,соответствующий индексам в массивах P_DNA и XY_DNA
                     n_lines_pod[j]++;//увеличение количества линий в j-том подсемействе
                     pod_exist=1;//подсемейство со значением 1-го параметра текущей (i-той) линии семейства D существует
                  }
               }
               if(pod_exist==0)//если подсемейства с таким значением первого параметра еще нет
               {
                  par_podsem[n_pod]=P_DNA[i_d][i][1];//значение параметра второго типа (с индексом 1) n_pod-го подсемейства семейства D
            										 //(равен параметру второго типа i-й линии семейства D)
               	index_lin_podsem[n_pod][0]=i;//i - индекс 0-й линии нового (n_pod-го) подсемейства,соответствующий индексам в массивах P_DNA и XY_DNA
                  n_lines_pod[n_pod]=1;//количества линий в новом подсемействе - 1
                  n_pod++;//увеличение количества подсемейств
               }
	         }
 	      	par_sem=(double**)calloc(n_pod,sizeof(double*));//указатель на массив значений параметров внутри подсемейства
	         for(i=0;i<=n_pod-1;i++)//перебрать все подсемейства
	         {
	            par_sem[i]=(double*)calloc(n_lines_pod[i],sizeof(double));//массив значений параметров внутри подсемейства
	         	for(j=0;j<=n_lines_pod[i]-1;j++)//перебрать все линии i-го подсемейства
 	           	{
	            	par_sem[i][j]=P_DNA[i_d][index_lin_podsem[i][j]][0];
              	}
	         }
            for(k=0; k<18; k++)//по всем рабочим длинам волн
            {
            	e=(double*)calloc(n_pod,sizeof(double));//указатель на массив значений D для каждого подсемейства
	         	for(i=0;i<=n_pod-1;i++)//по всем подсемействам
	         	{
	         		e1=(double*)calloc(n_lines_pod[i],sizeof(double));//указатель на массив значений D для каждой линии подсемейства
	            	for(j=0;j<=n_lines_pod[i]-1;j++)//по всем линиям подсемейства
                     //Интерполяция значений D для каждой линии внутри i-го подсемейства по значениям абсцисс для каждой рабочей длины волны
	            		e1[j]=linterp(DNA[Index_N_tochek[i_d][index_lin_podsem[i][j]]],XY_DNA[i_d][index_lin_podsem[i][j]][0],XY_DNA[i_d][index_lin_podsem[i][j]][1],Abscissa[k][DNA[Index_N_parametrov[i_d]-2]]);
                  //Интерполяция значений D внутри i-го подсемейства по первому (с индексом 0) параметру
	            	e[i]=linterp(n_lines_pod[i],par_sem[i],e1,Parameter[k][DNA[Index_N_parametrov[i_d]+1]]);
	            	free(e1);
	         	}
               //Интерполяция значений D между подсемействами по второму (с индексом 1) параметру
	         	D_m[k]=linterp(n_pod,par_podsem,e,Parameter[k][DNA[Index_N_parametrov[i_d]+2]]);
	         	free(e);
            }
	         for(i=0;i<=n_pod-1;i++)
	         	free(par_sem[i]);
	         free(par_sem);
	         /*free(nomer_first_lin_sem);*/
	         free(par_podsem);
	         free(n_lines_pod);
            for(i=0; i<=DNA[Index_N_liniy[i_d]]-1;i++)//освободить память по количеству линий семейства D начиная
            	free(index_lin_podsem[i]);
            free(index_lin_podsem);
         }
	      else if(DNA[Index_N_parametrov[i_d]]==1)//если количество параметров в семействе D - 1
	      {
	         par_sem=(double**)calloc(1,sizeof(double*));
	        	par_sem[0]=(double*)calloc(DNA[Index_N_liniy[i_d]],sizeof(double));
	         for(i=0;i<=DNA[Index_N_liniy[i_d]]-1;i++)
	         	par_sem[0][i]=P_DNA[i_d][i][0];//массив значений 0-го параметря для i-тых линий семейства D
            for(k=0; k<18; k++)//по всем рабочим длинам волн
            {
	         	e=(double*)calloc(DNA[Index_N_liniy[i_d]],sizeof(double));
	         	for(i=0;i<=DNA[Index_N_liniy[i_d]]-1;i++)
	         		e[i]=linterp(DNA[Index_N_tochek[i_d][i]],XY_DNA[i_d][i][0],XY_DNA[i_d][i][1],Abscissa[k][DNA[Index_N_parametrov[i_d]-2]]);
	         	D_m[k]=linterp(DNA[Index_N_liniy[i_d]],par_sem[0],e,Parameter[k][DNA[Index_N_parametrov[i_d]+1]]);
	         	free(e);
            }
	         free(par_sem[0]);
	         free(par_sem);
	      }
	      else//если количество параметров в семействе D - 0 (1 линия)
            for(k=0; k<18; k++)//по всем рабочим длинам волн
	      		D_m[k]=linterp(DNA[Index_N_tochek[i_d][0]],XY_DNA[i_d][0][0],XY_DNA[i_d][0][1],Abscissa[k][DNA[Index_N_parametrov[i_d]-2]]);
        }

        //Определение значения нормированной диаграммы направленности в горизонтальной плоскости (Е/Em)
        for(k=0; k<18; k++)//по всем рабочим длинам волн
        	e_em[k]=(double*)calloc(n_e_em,sizeof(double));//указатели на массивы с элементами по количеству семейств Е/Em
        for(o=0;o<=n_e_em-1;o++)
        {
        	if(DNA[Index_N_parametrov[i_e_em[o]]]>1)//если количество типов параметров в семействе Е/Em больше 1, т.е. 2
        	{
        		n_pod=1;//количество подсемейств (подсемейства группируются по равенству второго (с индексом 1) типа параметра семейства).
            		  //начальное значение - 1 подсемейство
        		n_lines_pod=(int*)calloc(DNA[Index_N_liniy[i_e_em[o]]],sizeof(int));//массив значений количества линий в подсемействе семейства Е/Em
        		n_lines_pod[0]=1;//количество линий в 0-м подсемействе семейства Е/Em (начальное значение - 1 линия в 0-м подсемействе)
        		par_podsem=(double*)calloc(DNA[Index_N_liniy[i_e_em[o]]],sizeof(double));//массив значений параметров второго типа (с индексом 1)
            																								 //подсемейства семейства Е/Em
        		par_podsem[0]=P_DNA[i_e_em[o]][0][1];//значение параметра второго типа (с индексом 1) 0-го подсемейства семейства Е/Em
            										 		 //(равен параметру второго типа 0-й линии семейства Е/Em)
            index_lin_podsem=(int**)calloc(DNA[Index_N_liniy[i_e_em[o]]],sizeof(int*));//массив индексов линий подсемейств семейства Е/Em,
            																							//соответствующих индексам в массивах P_DNA и XY_DNA
            for(i=0; i<=DNA[Index_N_liniy[i_e_em[o]]]-1;i++)//выделить память по количеству линий семейства Е/Em
            	index_lin_podsem[i]=(int*)calloc(DNA[Index_N_liniy[i_e_em[o]]],sizeof(int));//массив индексов линий i-го подсемейства семейства Е/Em,
            																									 //соответствующих индексам линий в массивах P_DNA и XY_DNA
            index_lin_podsem[0][0]=0;//индекс 0-ой линий (второй индекс) 0-го подсемейства (первый индекс) семейства Е/Em,
            								 //соответствующий индексам линий в массивах P_DNA и XY_DNA
         	for(i=1;i<=DNA[Index_N_liniy[i_e_em[o]]]-1;i++)//перебрать линии семейства Е/Em начиная с 1-й (информация о 0-й линии уже занесена в массивы)
         	{
            	pod_exist=0;//признак, что подсемейство со значением 1-го параметра текущей (i-той) линии семейства Е/Em существует
               				//0 - не существует, 1 - существует
               for(j=0; j<n_pod; j++)//Перебрать все найденные подсемейства семейства Е/Em
               {
               	if(P_DNA[i_e_em[o]][i][1]==P_DNA[i_e_em[o]][index_lin_podsem[j][0]][1])//Если значение 1-го параметра i-ой линии семейства Е/Em == значению
                  																						//1-го параметра 0-ой линии j-го подсемейства (т.е. подсемейство
                                                                                    //с таким значением 1-го параметра уже есть)
                  {
                     index_lin_podsem[j][n_lines_pod[j]]=i;//индекс текущей линии j-го подсемейства,соответствующий индексам в массивах P_DNA и XY_DNA
                     n_lines_pod[j]++;//увеличение количества линий в j-том подсемействе
                     pod_exist=1;//подсемейство со значением 1-го параметра текущей (i-той) линии семейства Е/Em существует
                  }
               }
               if(pod_exist==0)//если подсемейства с таким значением первого параметра еще нет
               {
                  par_podsem[n_pod]=P_DNA[i_e_em[o]][i][1];//значение параметра второго типа (с индексом 1) n_pod-го подсемейства семейства Е/Em
            										 //(равен параметру второго типа i-й линии семейства Е/Em)
               	index_lin_podsem[n_pod][0]=i;//i - индекс 0-й линии нового (n_pod-го) подсемейства,соответствующий индексам в массивах P_DNA и XY_DNA
                  n_lines_pod[n_pod]=1;//количества линий в новом подсемействе - 1
                  n_pod++;//увеличение количества подсемейств
               }
         	}
         	par_sem=(double**)calloc(n_pod,sizeof(double*));//указатель на массив значений параметров 0-го типа внутри подсемейства
         	for(i=0;i<=n_pod-1;i++)//перебрать все подсемейства
         	{
            	par_sem[i]=(double*)calloc(n_lines_pod[i],sizeof(double));//массив значений параметров 0-го типа внутри подсемейства
         		for(j=0;j<=n_lines_pod[i]-1;j++)//перебрать все линии i-го подсемейства
            	{
            		par_sem[i][j]=P_DNA[i_e_em[o]][index_lin_podsem[i][j]][0];
            	}
         	}
            for(k=0; k<18; k++)//по всем рабочим длинам волн
            {
         		e=(double*)calloc(n_pod,sizeof(double));//указатель на массив значений Е/Em для каждого подсемейства
         		for(i=0;i<=n_pod-1;i++)//по всем подсемействам
         		{
         			e1=(double*)calloc(n_lines_pod[i],sizeof(double));//указатель на массив значений Е/Em для каждой линии подсемейства
            		for(j=0;j<=n_lines_pod[i]-1;j++)//по всем линиям подсемейства
                     //Интерполяция значений Е/Em для каждой линии внутри i-го подсемейства по значениям абсцисс для каждой рабочей длины волны
            			e1[j]=linterp(DNA[Index_N_tochek[i_e_em[o]][index_lin_podsem[i][j]]],XY_DNA[i_e_em[o]][index_lin_podsem[i][j]][0],XY_DNA[i_e_em[o]][index_lin_podsem[i][j]][1],Abscissa[k][DNA[Index_N_parametrov[i_e_em[o]]-2]]);
                  //Интерполяция значений Е/Em внутри i-го подсемейства по первому (с индексом 0) параметру
            		e[i]=linterp(n_lines_pod[i],par_sem[i],e1,Parameter[k][DNA[Index_N_parametrov[i_e_em[o]]+1]]);
            		free(e1);
         		}
               //Интерполяция значений Е/Em между подсемействами по второму (с индексом 1) параметру для каждой рабочей длины волны
         		e_em[k][o]=linterp(n_pod,par_podsem,e,Parameter[k][DNA[Index_N_parametrov[i_e_em[o]]+2]]);
         		free(e);
            }
         	for(i=0;i<=n_pod-1;i++)
         		free(par_sem[i]);
         	free(par_sem);
         	free(par_podsem);
         	free(n_lines_pod);
            for(i=0; i<=DNA[Index_N_liniy[i_e_em[o]]]-1;i++)//освободить память по количеству линий семейства Е/Em начиная
            	free(index_lin_podsem[i]);
            free(index_lin_podsem);
        	}
        	else if(DNA[Index_N_parametrov[i_e_em[o]]]==1)//если количество параметров в семействе Е/Em - 1
        	{
         	par_sem=(double**)calloc(1,sizeof(double*));
        		par_sem[0]=(double*)calloc(DNA[Index_N_liniy[i_e_em[o]]],sizeof(double));
         	for(i=0;i<=DNA[Index_N_liniy[i_e_em[o]]]-1;i++)
         		par_sem[0][i]=P_DNA[i_e_em[o]][i][0];//массив значений 0-го параметря для i-тых линий семейства Е/Em
            for(k=0; k<18; k++)//по всем рабочим длинам волн
            {
         		e=(double*)calloc(DNA[Index_N_liniy[i_e_em[o]]],sizeof(double));
         		for(i=0;i<=DNA[Index_N_liniy[i_e_em[o]]]-1;i++)
         			e[i]=linterp(DNA[Index_N_tochek[i_e_em[o]][i]],XY_DNA[i_e_em[o]][i][0],XY_DNA[i_e_em[o]][i][1],Abscissa[k][DNA[Index_N_parametrov[i_e_em[o]]-2]]);
         		e_em[k][o]=linterp(DNA[Index_N_liniy[i_e_em[o]]],par_sem[0],e,Parameter[k][DNA[Index_N_parametrov[i_e_em[o]]+1]]);
         		free(e);
            }
         	free(par_sem[0]);
         	free(par_sem);
        	}
        	else//если количество параметров в семействе Е/Em - 0 (1 линия)
         	for(k=0; k<18; k++)//по всем рабочим длинам волн
        			e_em[k][o]=linterp(DNA[Index_N_tochek[i_e_em[o]][0]],XY_DNA[i_e_em[o]][0][0],XY_DNA[i_e_em[o]][0][1],Abscissa[k][DNA[Index_N_parametrov[i_e_em[o]]-2]]);
        }
        if(n_e_em!=0)
        	for(k=0; k<18; k++)//по всем рабочим длинам волн
        		for(i=0;i<=n_e_em-1;i++)
         		E_Em[k]=E_Em[k]*e_em[k][i];

        for(k=0; k<18; k++)//по всем рабочим длинам волн
        	free(e_em[k]);//освобождение памяти под массивы с элементами по количеству семейств Е/Em

        if(Pd_Pm_int==0)//если передающая антенна - коэффициент усиления
         for(k=0; k<18; k++)//по всем рабочим длинам волн
        		E_D->hReal[0][k]=(E_m[k]*E_Em[k]>0.0000316)?E_m[k]*E_Em[k]:0.0000316;//С учетом максимального ослабления в провалах ДНА 45 дБ(0.0000316 по напряженности)
        else//если приемная антенна - КНД
         for(k=0; k<18; k++)//по всем рабочим длинам волн
        		E_D->hReal[0][k]=(D_m[k]*E_Em[k]>0.0000316)?D_m[k]*E_Em[k]:0.0000316;//С учетом максимального ослабления в провалах ДНА 45 дБ(0.0000316 по напряженности)

        //Тип антенны по буквам
        j=strlen(DNA_KVA_name);//количество символов в типе антенны
        for(i=0; i<=j; i++ )
        {
        	E_D->hReal[1][i]=(double)((DNA_KVA_name[i]>=0)?DNA_KVA_name[i]:DNA_KVA_name[i]+256);//256 - смещение для русских букв
        }


        return 0;               // return 0 to indicate there was no error

    }// konec FDNAFunction





    char *FReadDNA_KVAFunction(double *DNA_KVA,
    										int Pd_Pm)			 // Читает файл диаграммы направленности коротковолновой антенны и построчно копирует его в массив
                                 						 //DNA_KVA.
                                                    // Pd_Pm - 0 - передающая антенна, 1 - приемная антенна.
                                                    // Возвращает указатель на строку, содержащую тип антенны
                                                    //В случае ошибки возвращает NULL

    {
      OPENFILENAME ofn;       // common dialog box structure
		TCHAR szFile[_MAX_PATH]; // buffer for filename      char DNA_KVA_name[256], s[256];//указатель на строку для чтения символных строк из файла		HWND hwnd;              // owner window		LPTSTR lpProcessCmdLine;      static char process_dir_DNA_KVA_file[_MAX_PATH];
      FILE *cur_dir_file, *DNA_KVA_stream;
      char DNA_KVA_PathFile_Name[_MAX_PATH];//Путь и имя выбранного файла диаграммы направленности
      double DNA_KVA_double[4096];
      int i, k, n_str=2;//счетчик считанных строк из файла диаграммы направленности (2-учитывает первые 2 строки - "0" и тип антенны)
      DWORD ofn_error;

      //считываем директорию исполняемого файла нашего процесса (Mathcad)
      lpProcessCmdLine=GetCommandLine();
      k=lstrlen(lpProcessCmdLine);//количество символов в строке пути включая имя исполняемого файла
      while(lpProcessCmdLine[k]!='\\') k--;//доходим до последнего слэша в пути\имени файла
      if(lstrcpyn(process_dir_DNA_KVA_file,(lpProcessCmdLine[0]=='\"')? &lpProcessCmdLine[1]:lpProcessCmdLine, k++))//копируем путь директории в буфер
      {
      	if(access(strcat(process_dir_DNA_KVA_file,"\\cur_DNA_KVA_dir"), 0)==0)//если файл cur_DNA_KVA_dir существует
         {
         	cur_dir_file=fopen(process_dir_DNA_KVA_file,"rt");
            fread(DNA_KVA_PathFile_Name,_MAX_PATH+1,1,cur_dir_file);
            fclose(cur_dir_file);
         }
         else
         	DNA_KVA_PathFile_Name[0]='\0';//если файла cur_DNA_KVA_dir не существует - обнуляем начальный путь
      }


      ZeroMemory(&ofn, sizeof(ofn));
    	ZeroMemory(szFile, sizeof(TCHAR)*_MAX_PATH);

     	// Initialize OPENFILENAME
     	ofn.lStructSize = sizeof(ofn);
     	hwnd=GetForegroundWindow();//хэндл текущего окна (Mathcad)
     	ofn.hwndOwner = hwnd;
     	ofn.lpstrFile = szFile;
     	ofn.nMaxFile = _MAX_PATH;
     	ofn.lpstrFilter = "Диаграммы направленности\0*.kva\0All\0*.*\0Text\0*.TXT\0";
     	ofn.nFilterIndex = 1;
     	ofn.lpstrFileTitle = NULL;
     	ofn.nMaxFileTitle = 0;
     	ofn.lpstrInitialDir = DNA_KVA_PathFile_Name;
      ofn.lpstrTitle = (Pd_Pm==0)?"Передающая антенна":"Приемная антенна";
     	ofn.Flags = OFN_EXPLORER;
      // Display the Open dialog box.		if (GetOpenFileName(&ofn)==TRUE)    	{      	cur_dir_file=fopen(process_dir_DNA_KVA_file,"wt");
         fwrite(szFile,strlen(szFile)+1,1,cur_dir_file);//запись пути директории, в которой последний раз открывалась диаграмма направленности КВ антенн,
         															  //в директорию Mathcad
         fclose(cur_dir_file);
      }      else      {      	//ofn_error=CommDlgExtendedError();      	return NULL;      }
      //Открываем поток для чтения файла диаграммы направленности
      DNA_KVA_stream=fopen(szFile,"rt");
      fgets(s,256,DNA_KVA_stream);//чтение 0-й строки из файла DNA_KVA ("0")
      fgets(DNA_KVA_name,256,DNA_KVA_stream);//чтение 1-й строки из файла DNA_KVA (тип антенны)
      DNA_KVA_name[strlen(DNA_KVA_name)-1]='\0';//удаление символа '\n'


      DNA_KVA[0]=0;//вместо 0-й строки из файла DNA_KVA
      DNA_KVA[1]=0;//вместо 1-й строки из файла DNA_KVA (тип антенны)

      while (!feof(DNA_KVA_stream))//пока не достигнут конец файла
      {
			fgets(s,256,DNA_KVA_stream);//чтение n_str-й строки из файла DNA_KVA ("0")
         DNA_KVA[n_str]=atof(s);//преобразование считанной строки в число и запись в массив
      	n_str++;//следующая строка
      }

      //Закрываем поток файла диаграммы направленности
      fclose(DNA_KVA_stream);

      return  &(DNA_KVA_name[0]);

    }// konec FReadDNA_KVAFunction



































/*
    FUNCTIONINFO    FAzRasst =
    {
    "FAzRasst",                          // Name by which mathcad will recognize the function
    "brad_a,lrad_a,brad_b,lrad_b",       // FAzRasst will be called as FAzRasst(brad_a,lrad_a,brad_b,lrad_b)
    "Vychislyaet azimut i rasstoyanie po koordinatam 2-h tochek",      // description of FAzRasst(brad_a,lrad_a,brad_b,lrad_b)
    (LPCFUNCTION)FAzRasstFunction,       // pointer to the executible code
    COMPLEX_ARRAY,                     // the return type is also a complex array
    4,                                  //  the function takes on 4 arguments
    { COMPLEX_SCALAR,COMPLEX_SCALAR,COMPLEX_SCALAR,COMPLEX_SCALAR}   // arguments are complex array and scalar
    };






























    LRESULT FAzRasstFunction(COMPLEXARRAY * const AzRasst,							// Opredelyaet azimut i rasstoyanie mejdu dvumya tochkami
                            COMPLEXSCALAR * const brad_a,						// AzRasst - 0-j element - pryamoj azimut v radianah, 1-j element - rasstojanie v km,
                            COMPLEXSCALAR * const lrad_a,
                            COMPLEXSCALAR * const brad_b,
                            COMPLEXSCALAR * const lrad_b) 						// brad_a - shirota v radianah tochki a,
                             																// lrad_a - dolgota v radianah tochki a,
                                                                           // brad_b - shirota v radianah tochki b,
                             																// lrad_b - dolgota v radianah tochki b

    {
		double eps=0.0000000000000005,e2=0.0066934216,z1=0.9966476701,w1,w2,sinu1,
      		 sinu2,cosu1,cosu2,L,x1,x2,v1,v2,dl1=0,dl2,ld,P,Q,a0,a1,sins,coss,s1,
             sina0,cos2a0,x,al,bt,az,bz,y,S;

      w1=sqrt(1-e2*sin(brad_a->real)*sin(brad_a->real));
      w2=sqrt(1-e2*sin(brad_b->real)*sin(brad_b->real));
      sinu1=sin(brad_a->real)*z1/w1;
      sinu2=sin(brad_b->real)*z1/w2;
      cosu1=cos(brad_a->real)/w1;
      cosu2=cos(brad_b->real)/w2;
      L=lrad_b->real-lrad_a->real;
      x1=sinu1*sinu2;
      x2=cosu1*cosu2;
      v1=cosu1*sinu2;
      v2=cosu2*sinu1;

      for(;;)
      {
			ld=L+dl1;
         P=cosu2*sin(ld);
         Q=v1-v2*cos(ld);
         if(P<0)
	         a1=2*M_PI+atan2(Q,P);
         else
            a1=atan2(Q,P);
         sins=M_PI*sin(a1)+Q*cos(a1);
         coss=x1+x2*cos(ld);
         s1=fabs(sins/coss);
         if(coss>0)
         	S=atan(s1);
         else
         	S=M_PI-atan(s1);

         sina0=cosu1*sin(a1);
         a0=asin(sina0);
         cos2a0=cos(a0)*cos(a0);
         x=2*x1-cos2a0*cos(S);
         al=(33523299-(28189-70*cos2a0)*cos2a0)*0.0000000001;
         bt=(28189-94*cos2a0)*0.0000000001;
         dl2=(al*S-bt*x*sin(S))*sina0;
         if(fabs(dl1-dl2)>eps)
         	dl1=dl2;
         else
         	break;
      }

      az=6356863.020+(10708.949-13.474*cos2a0)*cos2a0;
      bz=10708.938-17.956*cos2a0;
      y=(cos2a0*cos2a0-2*x*x)*cos(S);

      // allocate space for the return array Relyef
      if ( !MathcadArrayAllocate( AzRasst,  // allocate space for AzRasst
                    2,    //with 2 rows
                    1,    //and 2 cols
                    TRUE,   //  allocate the real part
                    FALSE    //  don't allocate the imaginary part
                    ) )
      return 2;           // if allocation is insufficient
                          // return the error code

      AzRasst->hReal[0][0]=a1;
      AzRasst->hReal[0][1]=az*S+(bz*x+4.487*y)*sin(S);


        return 0;               // return 0 to indicate there was no error

    }// konec FAzRasstFunction

*/



