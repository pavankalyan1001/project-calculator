//calculator
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<conio.h>
#define pi 3.14159265359
#define e 2.71828182846
int i, j,R,C;
int mat[10][10];
void quad();
void cubic();
void twroots();
void throots();
void comop();
void sumd();
void maxmin();
void median();
void matmul();
void mattra();
void echelon();
void inverse(int di);
float determinant(float [][25], float);
void cofactor(float [][25], float);
void transpose(float [][25], float [][25], float);
void display( int, int);
void input( int, int);
int Rank_Mat(int , int);
void swap(int, int, int);
void matrank();
int factorial(int n);
void quad()
{
	float a,b,c,x1,x2,dis,nume;
	printf("\n Enter the values of a,b,c in the form ax^2+bx+c=0:");
	scanf("%f%f%f",&a,&b,&c);
	dis=(b*b)-4*a*c; 
	if(dis<0)
	{
		dis=dis*-1;
		printf("\nThe roots of the qudratic equations:-%f+-%fi/%f",b,dis,2*a);
	}
	else if(dis>=0)
	{
	x1=(-b+sqrt(dis))/2*a;
	x2=(-b-sqrt(dis))/2*a;
	printf("\nThe roots of the quadratic equation:%f and %f",x1,x2);
	}
}
void cubic()
{
	double a,b,c,d,f,g,h,i,j,k,l,m,n,p,r,s,t,u,x1,x2,x3;
int w;
printf("\n Enter the values of a,b,c,d from the equation in the form a*x^3+b*x^2+c*x+d:\n");
scanf("%lf%lf%lf",&a,&b,&c);
f=((3*c/a)-(b*b/(a*a)))/3;
g=((2*b*b*b/(a*a*a))-(9*b*c/(a*a))+(27*d/a))/27; 
h=(g*g/4)+(f*f*f/27);
i=sqrt(((g*g/4)-h));
j=exp(log10(i)/log10(e)/3);
k=acos((-1)*(g/(2*i)));
l=j*(-1);
m=cos(k/3);
n=sqrt(3)*sin(k/3);
p=(b/3*a)*(-1);
r=(-1)*(g/2)+sqrt(h);
s=exp(log10(r)/log10(e)/3);
t=(-1)*(g/2)-sqrt(h);
u=exp(log10(t)/log10(e)/3);
if (h>0) w=1;
if (h<=0) w=3;
if ((f==0) && (g==0) && (h==0)) w=2;
switch (w){
case 1:
x1=(s+u)-(b/3*a);
x2=(-1)*(s+u)/2-(b/3*a);
x3=(s-u)*sqrt(3)/2;
printf("\nA 3 point:\n%lf\n%lf+i*%lf\n%lf -i*%lf", x1, x2, x3, x2, x3);
break;
case 2:
x1=exp(log10(d/a)/log10(e)/3)*(-1);
printf("\n There is a line:\n%lf", x1);
break;
case 3:
x1=2*j*cos(k/3)-(b/3*a);
x2=l*(m+n)+p;
x3=l*(m-n)+p;
printf("\nA 3 Roots:\n%lf\n%lf\n%lf", x1, x2, x3);
break;
}
}
void twroots()
{
	float a1, b1, c1, a2, b2, c2,x,y;
    printf("Enter the values for the first equation.");
    printf("\nEnter the values of a1,b1,c1 for the equation of the form a1x+b1y=c1:");
    scanf("%f%f%f",&a1,&b1,&c1);
    printf("Enter the values for the second equation.");
	printf("\nEnter the values of a2,b2,c2 for the equation of the form a2x+b2y=c2:");
    scanf("%f%f%f",&a2,&b2,&c2);
if ((a1 * b2)-(b1 * a2)== 0)
{	
    printf("The system has no solution.");
}
else{
    x = ((c1*b2) - (b1*c2))/((a1*b2)-(b1*a2));
    y = ((a1*c2) - (c1*a2)) / ((a1*b2) - (b1*a2));
    printf("\nThe solution of the given two equations:%f and %f",x,y);
	}	
}
void throots()
{
	float a1,b1,c1,d1,a2,b2,c2,d2,a3,D,b3,c3,d3,x,y,z;
printf("\nThe equations of the form a1x+b1y+c1z+d1=0\na2x+b2y+c2z+d2=0\na3x+b3y+c3z+d3=0\n ");
printf("\nEnter the values of coefficients of first equation(a1,b1,c1,d1):");
scanf("%f%f%f%f",&a1,&b1,&c1,&d1);
printf("\nEnter the values of coefficients of second equation(a2,b2,c2,d2):");
scanf("%f%f%f%f",&a2,&b2,&c2,&d2);
printf("\nEnter the values of coefficients of third equation(a3,b3,c3,d3):");
scanf("%f%f%f%f",&a3,&b3,&c3,&d3);
 D = (a1*b2*c3+b1*a3*c2+c1*a2*b3)-(a1*c2*b3+b1*a2*c3+c1*b2*a3);
 x = ((b1*c3*d2+c1*b2*d3+d1*c2*b3)-(b1*c2*d3+c1*b3*d2+d1*b2*c3))/D;
 y = ((a1*c2*d3+c1*a3*d2+d1*a2*c3)-(a1*c3*d2+c1*a2*d3+d1*c2*a3))/D;
 z = ((a1*b3*d2+b1*a2*d3+d1*b2*a3)-(a1*b2*d3+b1*a3*d2+d1*a2*b3))/D;
printf("The solutions to the above three equations are :");
printf("  x = %f\ny = %f\nz = %f",x,y,z);	
}
void comop()
{
	float a1,a2,b1,b2;
	printf("\nEnter the values of a1,b1 for first complex number of the form a1+ib1");
	scanf("%f%f",&a1,&b1);
	printf("\nEnter the values of a2,b2 for second complex number of the form a2+ib2");
	scanf("%f%f",&a2,&b2);
	printf("\nADDITION:");
	printf("\nSum of two complex numbers is:%f+i%f",a1+a2,b1+b2);
	printf("\nSUBTRACTION:");
	printf("\nDifference of two complex numbers is:%f+i%f",a1-a2,b1-b2);
	printf("\nMULTIPLICATION:");
	printf("\nMultiplication of two complex numbers is:%f+i%f",a1*a2-b1*b2,a1*b2+b1*a2);
}
void sumd()
{
	float a[100],s=0;
	int i,n;
	printf("\nEnter how many numbers to be added:");
	scanf("%d",&n);
	printf("\nEnter the numbers to be added:");
	for(i=1;i<=n;i++)
	{	
		printf("\nEnter the element%d:",i);
		scanf("%f",&a[i]);
		s=s+a[i];
	}
	printf("\nSum of numbers:%f",s);
}
void maxmin()
{
	int a[1000],i,n,min,max;
    printf("Enter the number of elements to be taken:");
    scanf("%d",&n);
    printf("Enter the elements: ");
    for(i=0;i<n;i++)
    {
    	printf("\nEnter the element%d:",i+1);
        scanf("%d",&a[i]);
    }
    min=max=a[0];
    for(i=1;i<n;i++)
    {
        if(min>a[i])
			min=a[i];   
		if(max<a[i])
			max=a[i];       
    }
    printf("\nminimum of taken observations is : %d",min);
    printf("\nmaximum of taken observations is : %d",max);
}
void median()
{
	int i,j,n;
   float median,a[1000],t;
   printf("Enter the number of elements to be taken:\n");
   scanf("%d",&n);
   printf("Input %d values \n",n);
   for(i=1;i<=n;i++)
   {
   		printf("\nEnter the value of element%d:",i);
    	scanf("%f", &a[i]);
   }
   for(i=1;i<=n-1;i++)
   { 
      for(j=1;j<=n-i;j++) 
	  {
         if (a[j]<=a[j+1]) 
		 { 
            t = a[j];
            a[j] = a[j+1];
            a[j+1] = t;
         }
         else
         continue ;
      }
   }
   if ( n % 2 == 0)
    {
	   median = (a[n/2] + a[n/2+1])/2.0 ;
	}
   else
    {
   		median = a[n/2 + 1];
	}
   printf("\nobservations arranged in descending order:");
   for(i=1;i<=n;i++)
    {
     	printf("%f",a[i]);
	}
   printf("\n\nMedian of the taken observations is :%f\n", median);
}
void matmul()
{
	int a[5][5],b[5][5],mul[5][5],*p,r1,c1,r2,c2,i,j,*q,*s,k;
	printf("\nEnter the rows and columnsof the first matrix respectively:");
	scanf("%d%d",&r1,&c1);
	p=&a[0][0];
	printf("\nEnter the elements of first the matrix:");
	for(i=0;i<r1;i++)
	{
		for(j=0;j<c1;j++)
		{
			printf("\n Enter the element a%d%d:",i+1,j+1);
			scanf("%d",(p+(i*c1)+j));
		}
	}
	matm:
	printf("\nEnter the rows and columnsof the second matrix respectively:");
	scanf("%d%d",&r2,&c2);
	if(c1!=r2)
	{
		printf("\ncomlums of first matrix(c1) is not equal to rows of second matrix(r2)");
		goto matm;
	}
	printf("\nEnter the elements of second matrix:");
	q=&b[0][0];
	for(i=0;i<r2;i++)
	{
		for(j=0;j<c2;j++)
		{
			printf("\nEnter the element b%d%d:",i+1,j+1);
			scanf("%d",(q+(i*c2)+j));
		}
	}
	printf("\nMultiplication of two matrices:\n");
	s=&mul[0][0];
	for(i=0;i<r1;i++)
	{
		for(j=0;j<c2;j++)
		{
			*(s
			+(i*c1)+j)=0;
			for(k=0;k<c1;k++)
			{
				*(s+(i*c1)+j)=*(s+(i*c1)+j)+*(p+(i*c1)+k)*(*(q+(k*c1)+j));
			}
		}
	}
	for(i=0;i<r1;i++)
	{
		for(j=0;j<c2;j++)
		{
			printf("%d\t",*(s+(i*c1)+j));
		}
		printf("\n");
	}
}
void mattra()
{
	int a[5][5],*p,r,c,i,j;
	printf("\nEnter the size of the matrix(rows and columns) respectively:");
	scanf("%d%d",&r,&c);
	p=&a[0][0];
	printf("\nEnter the elements of the matrix:");
	for(i=0;i<r;i++)
	{
		for(j=0;j<c;j++)
		{
			printf("\nEnter the element%d%d:",i+1,j+1);
			scanf("%d",(p+(i*c)+j));
		}
	}
	printf("\nThe transpose of the matrix:\n");
	for(i=0;i<c;i++)
	{
		for(j=0;j<r;j++)
		{
			printf("%d\t",*(p+(j*r)+i));
		}
		printf("\n");
	}
}
void echelon()
{
	int row,col,i,j,k,a,b;
	float mat[80][80],tem;	
	printf("\nEnter the number of rows : ");
	scanf("%d",&row);
	printf("\nEnter the number of columns : ");
	scanf("%d",&col);
	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
		{
			printf("\nEnter the %d,%d element : ",i+1,j+1);
			scanf("%f",&mat[i][j]);
		}
	}
	printf("Your Matrix is :: \n");
	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
		{
			printf("%f\t",mat[i][j]);
		}
		printf("\n");
	}
	for(k=0;k<row;k++)
	{	
		if( (mat[k][k]) != 1)
		{
			float temp = mat[k][k];
			if(temp == 0)
				continue;
			for(j=0;j<col;j++)
			{
				mat[k][j] = ( (mat[k][j]) / temp );
			}
		}	
		for(i=k+1;i<row;i++)
		{
			tem = mat[i][k];
			for(j=k;j<col;j++)
				{
					
					mat[i][j] = mat[i][j] - ( mat[k][j] * tem );
				}
		}
		if(k==row-1)
			printf("Row Echelon form is : \n\n");
		else
			printf("Step %d\n\n",k+1);
		for(a=0;a<row;a++)
		{
			for(b=0;b<col;b++)
			{
				if(mat[a][b] == -0)
					mat[a][b] = 0;
				printf("%f\t",mat[a][b]);
			}
			printf("\n");
		}
	}
}
void inverse(int di)
{
	float a[25][25], k, d;
  int i, j;
  printf("Enter the order of the Matrix : ");
  scanf("%f", &k);
  printf("Enter the elements of %fX%f Matrix : \n", k, k);
  for (i = 0;i < k; i++)
    {
     for (j = 0;j < k; j++)
       {
       	printf("\nEnter the element%d%d:",i+1,j+1);
        scanf("%f", &a[i][j]);
        }
    }
  d = determinant(a, k);
  if(d<0)
  {
  	d=-1*d;
  }
  if(di==2)
  {
  	printf("\ndeterminent of the matrix is:%f",d);
  	goto dend;
  }
  if (d == 0)
   printf("\nInverse of Entered Matrix is not possible\n");
  else
   cofactor(a, k);
   dend:
   printf("\n");
}
float determinant(float a[25][25], float k)
{
  float s = 1, det = 0, b[25][25];
  int i, j, m, n, c;
  if (k == 1)
    {
     return (a[0][0]);
    }
  else
    {
     det = 0;
     for (c = 0; c < k; c++)
       {
        m = 0;
        n = 0;
        for (i = 0;i < k; i++)
          {
            for (j = 0 ;j < k; j++)
              {
                b[i][j] = 0;
                if (i != 0 && j != c)
                 {
                   b[m][n] = a[i][j];
                   if (n < (k - 2))
                    n++;
                   else
                    {
                     n = 0;
                     m++;
                     }
                   }
               }
             }
          det = det + s * (a[0][c] * determinant(b, k - 1));
          s = -1 * s;
          }
    }
    return (det);
}
void cofactor(float num[25][25], float f)
{
 float b[25][25], fac[25][25];
 int p, q, m, n, i, j;
 for (q = 0;q < f; q++)
 {
   for (p = 0;p < f; p++)
    {
     m = 0;
     n = 0;
     for (i = 0;i < f; i++)
     {
       for (j = 0;j < f; j++)
        {
          if (i != q && j != p)
          {
            b[m][n] = num[i][j];
            if (n < (f - 2))
             n++;
            else
             {
               n = 0;
               m++;
               }
            }
        }
      }
      fac[q][p] = pow(-1, q + p) * determinant(b, f - 1);
    }
  }
  transpose(num, fac, f);
}
void transpose(float num[25][25], float fac[25][25], float r)
{
  int i, j;
  float b[25][25], inverse[25][25], d;
  for (i = 0;i < r; i++)
    {
     for (j = 0;j < r; j++)
       {
         b[i][j] = fac[j][i];
        }
    }
  d = determinant(num, r);
  for (i = 0;i < r; i++)
    {
     for (j = 0;j < r; j++)
       {
        inverse[i][j] = b[i][j] / d;
        }
    }
   printf("\n\n\nThe inverse of matrix is : \n");
   for (i = 0;i < r; i++)
    {
     for (j = 0;j < r; j++)
       {
         printf("\t%f", inverse[i][j]);
        }
    printf("\n");
     }
}
void swap( int row1,int row2, int col)
{
    for( i = 0; i < col; i++)
    {
        int temp = mat[row1][i];
        mat[row1][i] = mat[row2][i];
        mat[row2][i] = temp;
    }
}
int Rank_Mat(int row1, int col1)
{
    int r, c;
    for(r = 0; r< col1; r++)
    {
        display(R,C);
        if( mat[r][r] )  
        for(c = 0; c < row1; c++)
            if(c != r)
            {
                float ratio = mat[c][r]/ mat[r][r];
                for( i = 0; i < col1; i++)
                    mat[c][i] -= ratio * mat[r][i];
            }
            else
                printf("\n");
        else
        {
            for(c =  r+1 ; c < row1;  c++)
                if (mat[c][r])
                {
                    swap(r,c,col1);
                    break ;
                }

            if(c == row1)
            {
                -- col1;

                for(c = 0; c < row1; c ++)
                    mat[c][r] = mat[c][col1];
            }
            --r;
        }
    }
    return col1;
}
void display( int row, int col)
{
    for(i = 0; i < row; i++)
    {
        for(j = 0; j < col; j++)
        {
            printf("  %d", mat[i][j]);
        }
        printf("\n");
    }
}
void input( int row, int col)
{
    int value;
    printf("\nEnter the elements into the matrix");
    for(i = 0 ; i< row; i++)
    {
        for(j = 0 ;  j<col; j++)
        {
            printf("Input Value for: %d: %d: ", i+1, j+1);
            scanf("%d",  &value);
            mat[i][j] = value;
        }
    }
}
void matrank()
{
	int rank;
    printf("\n Enter number of rows:");
    scanf("%d", &R);
    printf("\n Enter number of columns:");
    scanf("%d", &C);
    input(R, C);
    rank = Rank_Mat(R, C);
    printf("\n Rank of above matrix is : %d", rank);
}
int factorial(int n)
{
	int i,fact=1;        
    for(i=1;i<=n;i++)
	{    
      fact=fact*i;    
  	}    
  return fact;
}
void main()
{
	int n,opt,A,B,C,D,E,F,G,H,I,J,K,des,z,stri,itri,htri;
	float as,bs,itriv,htriv;
	head:
	printf("\t \t \t SCIENTIFIC CALCULATOR");
	printf("\n1)BASIC CALCULATOR.\n2)TRIGONOMETRIC FUNCTIONS");
	printf("\n3)COMPLEX NUMBER OPERATIONS\n4)MATRIX");
	printf("\n5)EXPONENTIAL AND LOG OPERATIONS\n6)STATISTICS");
	printf("\n7)VECTOR\n8)SOLVING THE EQUATIONS");
	printf("\n9)CONVERTOR\n10)SOME ALGEBRAIC FUNCTIONS");
	printf("\n11)BASIC GEOMETRY\n12)IMPORTANT CONSTANTS");
	printf("\nEnter the operation to be performed from the menu given above:");
	scanf("%d",&n);
	switch (n)
	{
		case 1:
			A:
			printf("\n\t\tBASIC CALCULATOR:");
			sc:
			printf("\n1)ADDITION\n2)SUBTRACTION\n3)MULTIPLICATION\n4)DIVISION\n5)MODULUS/REMAINDER.");
			printf("\n6)SUM OF N NUMBERS.\n7)X^2\n8)SQRT(X)(SQUARE ROOT)\n9)X^3\n10)CBRT(X)(CUBE ROOT)");
			printf("\nEnter the operation to be performed from the menu given above:");
			scanf("%d",&opt);
			if(opt>=1&&opt<=5)
			{
				printf("\nEnter two numbers:");
				scanf("%f%f",&as,&bs);
				if(opt==1)
				{
					printf("\nSum of two numbers %f and %f is:%f",as,bs,as+bs);
				}
				else if(opt==2)
				{
					printf("\nDifference of two numbers %f and %f is:%f",as,bs,as-bs);
				}
				else if(opt==3)
				{
					printf("\nMultiplication of two numbers %f and %f is:%f",as,bs,as*bs);
				}	
				else if(opt==4)
				{
					printf("\nDivision  of two numbers %f and %f is:%f",as,bs,as/bs);
				}
				else if(opt==5)
				{
					printf("\n remainder when %f is divided with %f is:%d",as,bs,(int)as%(int)bs);
				}
			}
			else if(opt==6)
			{
				sumd();
			}
			else if(opt>=7&&opt<=10)
			{	
				printf("\nEnter the value of x:");
				scanf("%f",&as);
				if(opt==7)
				{
					printf("\nSquare of %f is:%f",as,pow(as,2));
				}
				else if(opt==8)
				{
					printf("\nSquare root of %f is:%f",as,sqrt(as));
				}
				else if(opt==9)
				{
					printf("\nCube of %f is:%f",as,pow(as,3));
				}
				else if(opt==10)
				{
					printf("\nCuberoot of %f is:%f",as,cbrt(as));
				}
			}
			else
			{
				printf("\nEnter a valid number.....\nchoose again\n");
				goto sc;
			}
			goto z;
			break;
		case 2:
			B:
			printf("\n\t\tTRIGONOMETRIC OPERATIONS:");
			alltrigo:
			printf("\n1)MAIN TRIGONOMETRIC FUNCTIONS\n2)INVERSE TRIGINOMETRIC FUNCTIONS");
			printf("\n3)HYPERBOLIC TRIGONOMETRIC FUNCTIONS");
			printf("\nEnter the operation to be performed from the menu given above:");
			scanf("%d",&opt);
			if(opt==1)
			{	
				printf("\n\t\tMAIN TRIGONOMETRIC FUNCTIONS");
				strig:
				printf("\n1)SINX");
				printf("\n2)COSX");
				printf("\n3)TANX");
				printf("\n4)COSECX");
				printf("\n5)SECX");
				printf("\n6)COTX");
				float ang;
				int dt;
				printf("\nEnter the trigonometric function:");
				scanf("%d",&stri);
				printf("\nEnter the angle in degrees:");
				scanf("%f",&ang);
				ang=(ang*pi)/180;
				dt=(2*ang)/pi;
				if(stri==1)
				{
					printf("\nsinx:%f",sin(ang));
				}
				else if(stri==2)
				{
					printf("\ncosx:%f",cos(ang));
				}
				else if(stri==3)
				{
					if(dt%2!=0)
					{
					printf("\nTanx value is Undefined for this angle");
					}
					else
					{
					printf("\ntanx=:%f",tan(ang));
					}
				}
				else if(stri==4)
				{
					if(dt%2==0)
					{
						printf("\nCosecx is  Undefined for this angle");
					}
					else
					{
						printf("\nCosecx=%f",1/sin(ang));
					}
				}
				else if(stri==5)
				{
					if(dt%2!=0)
					{
					printf("\nSecx value is Undefined for this angle");
					}
					else
					{
					printf("\nSecx=:%f",1/cos(ang));
					}
				}
				else if(stri==6)
				{
					if(dt%2==0)
					{
						printf("\nCotx is  Undefined for this angle");
					}
					else
					{
						printf("\nCotx=%f",1/tan(ang));
					}
				}
				else
				{
					printf("\nEnter a valid number....\nselect again");
					goto strig;
				}
			}
			else if(opt==2)
			{
				printf("\n\t\tINVERSE TRIGINOMETRIC FUNCTIONS"); 
				allitrigo:
				printf("\n1)arc(SINX)");
				printf("\n2)arc(COSX)");
				printf("\n3)arc(TANX)");
				printf("\n4)arc(COSECX)");
				printf("\n5)arc(SECX)");
				printf("\n6)arc(COTX)");
				printf("\nEnter the inverse trigonometric function:");
				scanf("%d",&itri);
				if(itri==1||itri==2)
				{
					printf("\nEnter the value of x in the domain(-1,+1):");
					istrigo:
					scanf("%d",&itriv);
					if(itriv>1||itriv<-1)
					{
						printf("\nEnter the value which is in the domain(-1,+1)....");
						goto istrigo;
					}
					if(itri==1)
					{
						printf("\narc(sinx)=%f",asin(itriv)*180/pi);
					}
					if(itri==2)
					{
						printf("\narc(cosx)=%f",acos(itriv)*180/pi);
					}
				}
				else if(itri==3||itri==6)
				{
					printf("\nEnter the value of x:");
					scanf("%f",&itriv);
					if(itri==3)
					{
						printf("\narc(tanx)=%f",atan(itriv)*180/pi);
					}
					if(itri==6)
					{
						if(itriv>0){
						printf("\narc(cotx)=%f",atan(1/itriv)*180/pi);
						}
						else if(itriv<0)
						{
							printf("\narc(cotx)=%f",(atan(1/itriv)+pi)*180/pi);
						}
						else
						{
							printf("\narc(cotx)=%f",pi/2);
						}
					}
				}
				else if(itri==4||itri==5)
				{
					printf("\nEnter the value of x(except(-1,+1)):");
					ictrigo:
					scanf("%f",&itriv);
					if(itriv>-1&&itriv<1)
					{
						printf("\nEnter the value of x not in (-1,+1).. ");
						goto ictrigo;
					}
					if(itri==4)
					{
						printf("\narc(cosecx)=%f",asin(1/itriv)*180/pi);
					}
					if(itri==5)
					{
						printf("\narc(secx)=%f",acos(1/itriv)*180/pi);
					}
				}
				else
				{
					printf("\nEnter the valid number....\nchoose again");
					goto allitrigo;
				}				
			}
			else if(opt==3)
			{
				printf("\n\t\tHYPERBOLIC TRIGONOMETRIC FUNCTIONS");
				htrigo:
				printf("\n1)sinhx\n2)coshx\n3)tanhx\n4)cosechx\n5)sechx\n6)cothx\n");
				printf("\nEnter the choice from above:");
				scanf("%d",&htri);
				if(htri==1||htri==2||htri==3||htri==5)
				{
					printf("\nEnter the value of x:");
					scanf("%d",&htriv);
					if(htri==1)
					{
						printf("\nsinhx=%f",sinh(htriv));
					}
					if(htri==2)
					{
						printf("\ncoshx=%f",cosh(htriv));
					}
					if(htri==3)
					{
						printf("\ntanhx=%f",tanh(htriv));
					}
					if(htri==5)
					{
						printf("\nsinhx=%f",1/cosh(htriv));
					}
				}
				if(htri==4||htri==6)
				{
					chtri:
					printf("\nEnter the value of x (other than 0):");
					scanf("%d",&htriv);
					if(htriv==0)
					{
						goto chtri;
					}
					if(htri==4)
					{
						printf("\ncosechx=%f",1/sinh(htriv));
					}
					if(htri==6)
					{
						printf("\ncothx=%f",1/tanh(htriv));
					}
				}
				else
				{
					printf("\nEnter a valid number.....\nchoose again");
					goto htrigo;
				}
			}
			else
			{
				printf("\nEnter a valid number.....\nchoose again....");
				goto alltrigo;
			}
			goto z;
			break;
		case 3:
			C:
			printf("\n\t\tCOMPLEX NUMBERS:");
			cmplx:
			printf("\n1)ARGUMENT OF THE COMPLEX NUMBERS.\n2)CONJUGATE OF THE COMPLEX NUMBERS.");
			printf("\n3)REAL PART\n4)IMAGINARY PART\n5)ABSOLUTE/MODULUS OF THE COMPLEX NUMBER.\n6)ADD/SUB/MUL");
			printf("\nEnter the operation to be performed from the menu given above:");
			scanf("%d",&opt);
			if(opt==1)
			{
				printf("\n\t\tARGUMENT OF THE COMPLEX NUMBERS");
			}
			else if(opt==2)
			{
				printf("\n\t\tCONJUGATE OF THE COMPLEX NUMBERS");
			}
			else if(opt==3)
			{
				printf("\n\t\tREAL PART");
			}
			else if(opt==4)
			{
				printf("\n\t\tIMAGINARY PART");
			}
			else if(opt==5)
			{
				printf("\n\t\tABSOLUTE/MODULUS OF THE COMPLEX NUMBER");
			}
			else if(opt==6)
			{
				printf("\n\t\t ADD/SUB/MUL");
				comop();
				goto z;
			}
			float a,b,s,arg;
			printf("\nEnter the values of a,b for the complex number of the format a+ib:");
			scanf("%f%f",&a,&b);
			if(opt==1)
			{
				arg=atan(b/a);
				printf("\nThe argument of the complex number %f+i%f is(in radians):%f",a,b,arg);
				printf("\nThe argument of the complex number %f+i%f is(in degrees):%f",a,b,arg*57.29577951308);
			}
			else if(opt==2)
			{	
				if(b>=0)
				printf("\nThe conjugate of the complex number %f+i%f is:%f-i%f",a,b,a,b);
				else
				printf("\nThe conjugate of the complex number %f-i%f is:%f+i%f",a,b*-1,a,b*-1);
			}
			else if(opt==3)
			{
				printf("\nThe real part of the complex number %f+i%f is:%f",a,b,a);
			}
			else if(opt==4)
			{
				printf("\nThe imaginary part of the complex number %f+i%f is:%f",a,b,b);
			}
			else if(opt==5)
			{
				s=sqrt(pow(a,2)+pow(b,2));
				printf("\nThe absolute/mod value of the complex number %f+i%f is:%f",a,b,s);
			}
			if(n>6||n<=0)
			{
				printf("\nEnter a valid number....\nchoose again...");
				goto cmplx;
			}
			break;
			goto z;
		case 4:
			D:
			printf("\nMATRIX OPERATIONS:");
			mat:
			printf("\n1)ADDITION\n2)SUBTRACTION\n3)MULTIPLICATION\n4)TRANSPOSE OF THE MATRIX\n5)INVERSE OF THE MATRIX");
			printf("\n6)RANK OF MATRIX\n7)DETERMINANT OF THE MATRIX\n8)ECHELON FORM OF THE MATRIX");
			printf("\nEnter the operation to be performed from the menu given above:");
			scanf("%d",&opt);
			if(opt==1||opt==2)
			{
				int r, c, a[100][100], b[100][100],sum[100][100],i,j;
 				printf("Enter the number of rows: ");
  				scanf("%d", &r);
  				printf("Enter the number of columns: ");
  				scanf("%d", &c);
  				printf("\nEnter elements of 1st matrix:\n");
  				for (i = 0; i < r; ++i)
				{
    				for (j = 0; j < c; ++j)
					{
      					printf("Enter element a%d%d: ", i + 1, j + 1);
      					scanf("%d", &a[i][j]);
					}
				}
  				printf("Enter elements of 2nd matrix:\n");
  				for (i = 0; i < r; ++i){
    			for (j = 0; j < c; ++j) {
      			printf("Enter element b%d%d: ", i + 1, j + 1);
      			scanf("%d", &b[i][j]);
      		}
    		}
    		if(opt==1)
    		{
    			printf("\nADDITION OF MATRIX");
    			for (i = 0; i < r; ++i)
  				{
    			for (j = 0; j < c; ++j) 
				{
      				sum[i][j] = a[i][j] + b[i][j];
    			}	
				}
  				printf("\nSum of two matrices: \n");
  			for (i = 0; i < r; ++i)
  			{	
   			for (j = 0; j < c; ++j) 
				{
      		printf("%d\t", sum[i][j]); 
      			}
      		printf("\n");
			}
			}
			if(opt==2)
			{
				printf("\nSUBTRACTION OF MATRIX");
				for (i = 0; i < r; ++i)
  			{
    		for (j = 0; j < c; ++j) 
				{
      				sum[i][j] = a[i][j] - b[i][j];
    			}	
			}
  		printf("\nDifference of two matrices: \n");
  		for (i = 0; i < r; ++i)
  			{
    		for (j = 0; j < c; ++j) 
				{
      				printf("%d\t", sum[i][j]); 
      			}
      			printf("\n");
			}
			}
		}
			else if(opt==3)
			{
				printf("\n\t\tMULTIPLICATION");
				matmul();
			}
			else if(opt==4)
			{
				printf("\n\t\tTRANSPOSE OF THE MATRIX");
				mattra();
			}
			else if(opt==5)
			{
				printf("\n\t\tINVERSE OF THE MATRIX");
				inverse(1);
			}
			else if(opt==6)
			{
				printf("\n\t\tRANK OF MATRIX");
				matrank();
			}
			else if(opt==7)
			{
				printf("\n\t\tDETERMINANT OF THE MATRIX");
				inverse(2);
			}
			else if(opt==8)
			{
				printf("\n\t\tECHELON FORM OF THE MATRIX");
				echelon();
			}
			else
			{
				printf("\nEnter a valid number.....choose again");
				goto mat;
			}
			goto z;
			break;
		case 5:
			E:
			printf("\nEXPONENTIAL AND LOG OPERATIONS:");
			exp:
			printf("\n1)X^Y\n2)e^X\n3)10^X\n4)log(X)\n5)log10(X)\n6)Ae^x\n7)Ae^x+B");
			printf("\n8)Alogx\n9)Alogx+B\n10)Ae^Bx+C");
			printf("\nEnter the operation to be performed from the menu given above:");
			scanf("%d",&opt);
			float xe,ye,ae,be,ce;
			if(opt==1)
			{	
				printf("\nEnter the values of x and y respectively:");
				scanf("%f%f",&xe,&ye);
				printf("\nThe value of %f^%f is:%f",xe,ye,pow(xe,ye));
			}
			else if(opt>=2&&opt<=5)
			{
				printf("\nEnter the value of x:");
				scanf("%f",&xe);
				if(opt==2)
				{
					printf("\nThe value of e^%f is:%f",xe,exp(xe));
				}
				else if(opt==3)
				{
					printf("\nThe value of 10^%f is:%f",xe,pow(10,xe));
				}
				else if(opt==4)
				{
					printf("\nThe value of log(%f) is:%f",xe,log(xe));
				}
				else if(opt==5)
				{
					printf("\nThe value of log10(%f) is:%f",xe,log10(xe));
				}
			}
			else if(opt==6||opt==8)
			{
				printf("\nEnter the value of a:");
				scanf("%f",&ae);
				printf("\nEnter the value of x:");
				scanf("%f",&xe);
				if(opt==6)
				{
					printf("\nThe value of %fe^%f is:%f",ae,xe,ae*exp(xe));
				}
				else
				{
					printf("\nThe value of %flog(%f) is:%f",ae,xe,ae*log(xe));
				}
			}
			else if(opt==7||opt==9)
			{
				printf("\nEnter the value of a:");
				scanf("%f",&ae);
				printf("\nEnter the value of b:");
				scanf("%f",&be);
				printf("\nEnter the value of x:");
				scanf("%f",&xe);
				if(opt==7)
				{
					printf("\nThe value of %fe^%f+%f is:%f",ae,xe,be,ae*exp(xe)+be);
				}
				else
				{
					printf("\nThe value of %flog(%f)+%f is:%f",ae,xe,be,ae*log(xe)+be);
				}
			}
			else if(opt==10)
			{
				printf("\nEnter the value of a:");
				scanf("%f",&ae);
				printf("\nEnter the value of b:");
				scanf("%f",&be);
				printf("\nEnter the value of c:");
				scanf("%f",&ce);
				printf("\nEnter the value of x:");
				scanf("%f",&xe);
				printf("\nThe value of %fe^%f*%f+%f is:%f",ae,be,xe,ce,ae*exp(be*xe)+ce);
			}
			else
			{
				printf("\nEnter a valid number.....\nchoose again");
				goto exp;
			}
			goto z;
			break;
		case 6:
			F:
			printf("\nSTATISTICS:");
			stats:
			printf("\n1)MEAN OF N VALUES.\n2)MEDIAN OF N VALUES.\n3)STANDARD DEVIATION.\n4)VARIANCE)");
			printf("\n5)MINIMUM AND MAXIMUM OF N NUMBERS.");
			printf("\nEnter the operation to be performed from the menu given above:");
			scanf("%d",&opt);
			float xs[2000];
    		int  is,ns;
    		float average, variance, std_deviation, sum = 0, sum1 = 0;
    		if(opt==1||opt==3||opt==4)
			{printf("Enter the number of elements to be taken: \n");
    		scanf("%d",&ns);
    		printf("Enter %d real numbers \n",ns);
    		for (is=0;is<ns;is++)
    		{
    			printf("\nEnter the value of %d element:",is+1);
        		scanf("%f", &xs[is]);
    		}
    		for(is=0;is<ns;is++)
    		{
        		sum=sum+xs[is];
    		}
    		average=sum/(float)ns;
    		for(is=0;is<ns;is++)
    		{
        		sum1=sum1+pow((xs[is] - average),2);
    		}
    		variance=sum1/(float)ns;
    		std_deviation=sqrt(variance);
    		}
			if(opt==1)
			{
				printf("\n\t\tMEAN OF N VALUES");
				printf("\nAverage/Mean of all elements = %f\n",average);
			}
			else if(opt==2)
			{
				printf("\n\t\tMEDIAN OF N VALUES\n");
				median();
			}
			else if(opt==3)
			{
				printf("\n\t\tSTANDARD DEVIATION");
				printf("Standard deviation = %f\n",std_deviation);
			}
			else if(opt==4)
			{
				printf("\n\t\tVARIANCE");
				printf("\nvariance of all elements = %f\n",variance);
			}
			else if(opt==5)
			{
				printf("\n\t\tMINIMUM AND MAXIMUM OF N NUMBERS\n");
				maxmin();
			}
			else
			{
				printf("\nEnter a valid number.....choose again");
				goto stats;
			}
			goto z;
			break;
		case 7:
			G:
			printf("\nVECTOR OPERATIONS");
			vector:
			printf("\n1)DOT PRODUCT\n2)ANGLE B/N TWO VECTORS\n3)CROSS PRODUCT\n4)PROJECTION(VECTOR 1,VECTOR 2)");
			printf("\nEnter the operation to be performed from the menu given above:");
			scanf("%d",&opt);
				if(opt==1)
				{
					printf("\n\t\tDOT PRODUCT");
				}
				if(opt==2)
				{
					printf("\n\t\tANGLE B/N TWO VECTORS");
				}
				if(opt==3)
				{
					printf("\n\t\tCROSS PRODUCT");
				}
				if(opt==4)
				{
					printf("\n\t\tPROJECTION(VECTOR 1,VECTOR 2)");
				}
				float a1,a2,a3,b1,b2,b3,dot,cross,ang,pro,x,y,z,d,div,co,mb,ma;
				printf("\n Enter the values for first vector:");
				printf("\n Enter the values of a1,a2,a3 for a1i+a2j+a3k respectively:");
				scanf("%f%f%f",&a1,&a2,&a3);
				printf("\n Enter the values for second vector:");
				printf("\n Enter the values of b1,b2,b3 for b1i+b2j+b3k respectively:");
				scanf("%f%f%f",&b1,&b2,&b3);
				if(opt==1)
				{
				dot=(a1*b1)+(a2*b2)+(a3*b3);
				printf("\n the dot product for the vectors :%f",dot);
				}
				else if(opt==2)
				{
					dot=(a1*b1)+(a2*b2)+(a3*b3);
					d=(sqrt(pow(a1,2)+pow(a2,2)+pow(a3,2)))*(sqrt(pow(b1,2)+pow(b2,2)+pow(b3,2)));
					div=dot/d;
					co=acos(div);
					printf("\n angle b/n two vectors(in radians):%f",co);
					printf("\n angle b/n two vectors(in degrees):%f",co*57.29577951308);
				}
				else if(opt==3)
				{
					x=(a2*b3)-(b2*a3);
					y=(a1*b3)-(b1*a3);
					y=b*-1;
					z=(a1*b2)-(a2*b1);
					printf("\nThe cross product of the vectors:%fi+%fj+%fk",x,y,z);
				}
				else if(opt==4)
				{
					dot=(a1*b1)+(a2*b2)+(a3*b3);
					mb=sqrt(pow(b1,2)+pow(b2,2)+pow(b3,2));
					ma=sqrt(pow(a1,2)+pow(a2,2)+pow(a3,2));
					printf("\nThe projection of a on b is :%f",dot/mb);
					printf("\nThe projection of b on a is :%f",dot/ma);
				}
				if(opt>=5||opt<=0)
				{
					printf("\nEnter a valid number.....\nchoose again.....");
					goto vector;
				}
				goto z;
			break;
		case 8:
			H:
			printf("\nSOLVING THE EQUATIONS:");
			soleq:
			printf("\n1)ROOTS OF THE EQUATION ax^2+bx+c=0.\n2)ROOTS OF THE EQUATION ax^3+bx^2+cx+d=0.");
			printf("\n3)ROOTS OF TWO EQUATIONS WITH TWO UNKNOWNS.\n4)ROOTS OF THREE EQUATIONS WITH THREE UNKNOWNS.");
			printf("\nEnter the operation to be performed from the menu given above:");
			scanf("%d",&opt);
			if(opt==1)
			{	printf("\n\t\tROOTS OF THE EQUATION ax^2+bx+c=0:");
				quad();
				goto z;
			}
			else if(opt==2)
			{	
				printf("\n\t\tROOTS OF THE EQUATION ax^3+bx^2+cx+d=0:");
				cubic();
				goto z;
			}
			else if(opt==3)
			{
				printf("\n\t\tROOTS OF TWO EQUATIONS WITH TWO UNKNOWNS:");
				twroots();
				goto z;
			}
			else if(opt==4)
			{
				printf("\n\t\tROOTS OF THREE EQUATIONS WITH THREE UNKNOWNS:");
				throots();
			}
			else
			{
				printf("\nEnter a valid number.....\nchoose again.....");
				goto soleq;
			}
			goto z;
			break;
		case 9:
			I:
			printf("---WELCOME TO CONVERSION CALCULATOR---\n");
			ocon:
			printf("Select what type of conversions you want to perform from the options below:\n");
			printf("1.Volume \t 2.Length \t 3.Weight and Mass \n");
			printf("4.Temperature \t 5.Energy \t 6.Area \n");
			printf("7.Speed \t 8.Power \t 9.Data \n");
			printf("10.Pressure \t 11.Angle \t 12.Bin,Oct,Dec,HexDec \n");
			printf("\nEnter the choice from the above list");
	int abcd;
	scanf("%d",&abcd);
		if(abcd==1){
			printf("---WELCOME TO VOLUME CONVERSIONS---\n");
			printf("Select the conversion:\n");
			ivcon:
			printf("1.litre -> gallon US \t 2.gallon US -> Litre \t 3.litre -> CC \t 4.CC -> Litre\n");
			int x;
			scanf("%d",&x);
			printf("Enter the operand:\n");
			float y;
			scanf("%f",&y);
			if(x==1)
			{
				float gal=0.264172*y;
				printf("%f gallons",gal);
			}
			if(x==2)
			{
				float lit=3.785412*y;
				printf("%f litres",lit);
			}
			if(x==3)
			{
				float cc=1000*y;
				printf("%f cubic centimeters",cc);
			}
			if(x==4)
			{
				float litr=0.001*y;
				printf("%f litres",litr);
			}
			else
			{
				printf("\nEnter a valid number.....\nchoose again.....");
				goto ivcon;
			}
		}
		else if(abcd==2)
		{
			printf("---WELCOME TO LENGTH CONVERSIONS---\n");
			printf("Select the conversion:\n");
			ilcon:
			printf("1. inch -> centimeter \t 2. centimeter -> inch \t 3. yard -> meter \t 4. meter -> yard\n");
			printf("5. feet -> meter \t 6. meter -> feet \t 7. mile -> kilometer \t 8. kilometer -> mile\n");
			int ab;
			scanf("%d",&ab);
			printf("Enter the operand:\n");
			float bc;
			scanf("%f",&bc);
			if(ab==1)
			{
				float cm=(bc)*2.54;
				printf("%f cm",cm);
			}
			if(ab==2)
			{
				float inc=(bc)/2.54;
				printf("%f inches",inc);
			}
			if(ab==3)
			{
				float mt=(bc)*0.9144;
				printf("%f meters",mt);
			}
			if(ab==4)
			{
				float yrd=(bc)/0.9144;
				printf("%f yards",yrd);
			}
			if(ab==5)
			{
				float met=(bc)*0.3048;
				printf("%f meters",met);
			}
			if(ab==6)
			{
				float fet=(bc)/0.3048;
				printf("%f feets",fet);
			}
			if(ab==7)
			{
				float kilm=(bc)/1.609344;
				printf("%f kilometers",kilm);
			}
			if(ab==8)
			{
				float mil=(bc)*2.54;
				printf("%f miles",mil);
			}
			else
			{
				printf("\nEnter a valid number.....\nchoose again.....");
				goto ilcon;
			}
	}
		else if(abcd==3){
		
			printf("---WELCOME TO WEIGHT AND MASS CONVERSIONS---\n");
			printf("Select the conversion:\n");
			iwmcon:
			printf("1. Oz -> grams \t 2. grams -> Oz \t 3. Lb -> Kg \t 4. Kg -> Lb\n");
			int xc;
			scanf("%d",&xc);
			printf("Enter the operand:\n");
			float zx;
			scanf("%f",&zx);
			if(xc==1)
			{
				float oz=(zx)*28.34952;
				printf("%f grams",oz);
			}
			if(xc==2)
			{
				float grms=(zx)/28.34952;
				printf("%f Oz",grms);
			}
			if(xc==3)
			{
				float kg=(zx)*0.453592;
				printf("%f kg",kg);
			}
			if(xc==4)
			{
				float lb=(zx)/0.453592;
				printf("%f Lb",lb);
			}
			else
			{
				printf("\nEnter a valid number.....\nchoose again.....");
				goto iwmcon;
			}
	}
		else if(abcd==4){
			printf("WELCOME TO TEMPERATURE CONVERSION\n");
			itemcon:
			printf("Convert The Following To Others \n 1.Celsius \t 2.Fahreinheit \t 3.Kelvin: \n");
			float bgt;
			scanf("%f",&bgt);
			float ert;
			printf("enter the operand:");
			scanf("%f",&ert);
			if(bgt==1)
			{
				printf("%f Celsius is:\n",ert);
				printf("%f Fahreinheit \n",((ert)*1.8)+32);
				printf("%f Kelvin \n",(ert)+273.15);
			}
			if(bgt==2)
			{
				printf("%f Fahreinheit is:\n",ert);
				printf("%f Celsius \n",(((ert)-32)*5)/9);
				printf("%f Kelvin \n",((((ert)-32)*5)/9)+273.15);
			}
			if(bgt==3)
			{
				printf("%f Kelvin is:\n",ert);
				printf("%f Celsius \n",(ert)-273.15);
				printf("%f Fahreinheit \n",((ert-273.15)*(1.8))+32);
			}
			else
			{
				printf("\nEnter a valid number.....\nchoose again.....");
				goto itemcon;
			}
		}
		else if(abcd==5){
			printf("---WELCOME TO ENERGY CONVERSIONS---\n");
			printf("Select the conversion:\n");
			iecon:
			printf("1. Joule -> calories \t 2. calories -> Joules \n");
			int qw;
			scanf("%d",&qw);
			printf("Enter the operand:\n");
			float we;
			scanf("%f",&we);
			if(qw==1)
			{
				float calr=(we)*0.239006;
				printf("%f Calories",calr);
			}
			if(qw==2)
			{
				float jol=(we)/0.239006;
				printf("%f Joule",jol);
			}
			else
			{
				printf("\nEnter a valid number.....\nchoose again.....");
				goto iecon;
			}
		}
		else if(abcd==6){
			printf("WELCOME TO AREA CONVERSION\n");
			iarcon:
			printf("1.Acre -> SquareFeet \t 2.SquareFeet -> Acre \t 3.Acre -> Hectare \t 4.Hectare -> Acre  \n ");
			printf("Select one from above to continue: ");
			float zxc;
			scanf("%f",&zxc);
			float qwe;
			printf("enter the operand:");
			scanf("%f",&qwe);
			if(zxc==1)
			{
				float sqft=(qwe)*0.000023;
				printf("%f SquareFeet",sqft);
			}
			if(zxc==2)
			{
				float acr=(qwe)/0.000023;
				printf("%f Acre",acr);
			}
			if(zxc==3)
			{
				float hect=(qwe)*0.404686;
				printf("%f Hectares",hect);
			}
			if(zxc==4)
			{
				float acre=(qwe)/0.404686;
				printf("%f Acre",acre);
			}
			else
			{
				printf("\nEnter a valid number.....\nchoose again.....");
				goto iarcon;
			}
		}
		else if(abcd==7){
			printf("WELCOME TO SPEED CONVERSION\n");
			ispcon:
			printf("1. Kilometers/hour -> Meters/Second \t 2.Meters/second -> KiloMeters/Hour \t 3.Miles/Hour -> Kilometers/Hour \t 4.Kilometers/Hour -> Miles/Hour  \n ");
			printf("Select one from above to continue: ");
			float rty;
			scanf("%f",&rty);
			float dfg;
			printf("enter the operand:");
			scanf("%f",&dfg);
			if(rty==1)
			{
				float mps=(dfg)*0.277778;
				printf("%f Meters/Second",mps);
			}
			if(rty==2)
			{
				float kmph=(dfg)/0.277778;
				printf("%f Kilometers/Hour",kmph);
			}
			if(rty==3)
			{
				float kkmph=(dfg)*1.6092;
				printf("%f Kilometers/Hour",kkmph);
			}
			if(rty==4)
			{
				float mph=(dfg)/1.6092;
				printf("%f Miles/Hour",mph);
			}
			else
			{
				printf("\nEnter a valid number.....\nchoose again.....");
				goto ispcon;
			}	
		}
		else if(abcd==8){
			printf("WELCOME TO POWER CONVERSION\n");
			ipocon:
			printf("1. KiloWatts -> Watts \t 2.Watts -> KiloWatts \t 3.HorsePower -> Watts \t 4.Watts -> HorsePower  \n ");
			printf("Select one from above to continue: ");
			float jkl;
			scanf("%f",&jkl);
			float cvb;
			printf("enter the operand:");
			scanf("%f",&cvb);
			if(jkl==1)
			{
				float watts=(cvb)*1000;
				printf("%f Watts",watts);
			}
			if(jkl==2)
			{
				float kwats=(cvb)/1000;
				printf("%f KiloWatts",kwats);
			}
			if(jkl==3)
			{
				float wats=(cvb)*745.6999;
				printf("%f Watts",wats);
			}
			if(jkl==4)
			{
				float hpwr=(cvb)/745.6999;
				printf("%f HorsePower",hpwr);
			}
			else
			{
				printf("\nEnter a valid number.....\nchoose again.....");
				goto ipocon;
			}
		}
		else if(abcd==9){
			printf("WELCOME TO DATA CONVERSION\n");
			idacon:
			printf("1. Megabytes -> GigaBytes \t 2.GigaBytes -> MegaBytes \t 3.MegaBits -> MegaBytes \t 4.MegaBytes -> Megabits \n ");
			printf("5. GigaBytes -> TerraBytes \t 6.TerraBytes -> GigaBytes \n ");
			printf("Select one from above to continue: ");
			float xcv;
			scanf("%f",&xcv);
			float ghj;
			printf("enter the operand:");
			scanf("%f",&ghj);
			if(xcv==1)
			{
				float gb=(ghj)*0.001;
				printf("%f GigaBytes",gb);
			}
			if(xcv==2)
			{
				float mb=(ghj)/0.001;
				printf("%f MegaBytes",mb);
			}
			if(xcv==3)
			{
				float mby=(ghj)*0.125;
				printf("%f MegaBytes",mby);
			}
			if(xcv==4)
			{
				float mbi=(ghj)/0.125;
				printf("%f MegaBits",mbi);
			}
			if(xcv==5)
			{
				float tb=(ghj)*0.001;
				printf("%f TeraBytes",tb);
			}
			if(xcv==6)
			{
				float gby=(ghj)/0.001;
				printf("%f GigaBytes",gby);
			}
			else
			{
				printf("\nEnter a valid number.....\nchoose again.....");
				goto idacon;
			}
		}
		else if(abcd==10){
			printf("WELCOME TO PRESSURE CONVERSION\n");
			iprcon:
			printf("Convert The Following To Others \n 1.Atmospheres \t 2.mmHg \t 3.Pascals \t 4.Bars : \n");
			float hjk;
			scanf("%f",&hjk);
			float qaz;
			printf("enter the operand:");
			scanf("%f",&qaz);
			if(hjk==1)
			{
				printf("%f Atmospheres is:\n",qaz);
				printf("%f mmHg \n",(qaz)*760.1275);
				printf("%f Pascals \n",(qaz)*101325);
				printf("%f Bars",(qaz)*1.01325);
			}
			if(hjk==2)
			{
				printf("%f mmHg is:\n",qaz);
				printf("%f Atmospheres \n",(qaz)*0.001316);
				printf("%f Pascals \n",(qaz)*133.3);
				printf("%f Bars",(qaz)*0.001333);
			}
			if(hjk==3)
			{
				printf("%f Pascals is:\n",qaz);
				printf("%f mmHg \n",(qaz)*0.007502);
				printf("%f Atmospheres \n",(qaz)*0.00001);
				printf("%f Bars",(qaz)*0.00001);
			}
			if(hjk==4)
			{
				printf("%f Bars is:\n",qaz);
				printf("%f mmHg \n",(qaz)*750.1875);
				printf("%f Pascals \n",(qaz)*100000);
				printf("%f Atmospheres",(qaz)*0.986923);
			}
			else
			{
				printf("\nEnter a valid number.....\nchoose again.....");
				goto iprcon;
			}
		}
		else if(abcd==11){
			printf("WELCOME TO ANGLE CONVERSION\n");
			iancon:
			printf("1. Degrees -> Radians \t 2.Radians -> Degrees  \n ");
			printf("Select one from above to continue: ");
			float yui;
			scanf("%f",&yui);
			float iop;
			printf("enter the operand:");
			scanf("%f",&iop);
			if(yui==1)
			{
				float radn=(iop)*0.017453;
				printf("%f Radians",radn);
			}
			if(yui==2)
			{
				float deg=(iop)*57.29578;
				printf("%f Degrees",deg);
			}
			else
			{
				printf("\nEnter a valid number.....\nchoose again.....");
				goto iancon;
			}
		}
		else if(abcd==12){
			printf("WELCOME TO BINARY CONVERSIONS\n");
			ibincon:
			printf("1. Binary -> Decimal \t 2.Decimal -> Binary \t 3.Binary -> HexaDecimal \t 4.HexaDecimal -> Binary \n ");
			printf("5. Binary -> Octal \t 6.Octal -> Binary   \n ");
			printf("Select one from above to continue: ");
			float bnm;
			scanf("%f",&bnm);
			if(bnm==1)
			{
				long long n;
				printf("Enter a binary number: ");
				scanf("%lld", &n);
				int dec = 0, i = 0, rem;
				while (n!=0) {
				rem = n % 10;
				 n /= 10;
				dec += rem * pow(2, i);
				++i;
				}
				printf("%lld in binary = %d in decimal", n,dec);		  
			}
			else if(bnm==2)
			{
				int na, bin;
  				printf("Enter a decimal number: ");
  				scanf("%d", &na);	
				long long bina = 0;
  				int rem, i = 1;
  				while (na!=0) {
    				rem = na % 2;
   					 na /= 2;
    				bina += rem * i;
   					 i *= 10;
  					}
  					printf("%lld in binary",bina);
	  
			}
			else if(bnm==3)
			{
				long int binaryval, hexadecimalval = 0, i = 1, remainder;
   				printf("Enter the binary number: ");
   				scanf("%ld", &binaryval);
   				while (binaryval != 0){
      				remainder = binaryval % 10;
      				hexadecimalval = hexadecimalval + remainder * i;
      				i = i * 2;
      				binaryval = binaryval / 10;
   						}
  				 printf("Equivalent hexadecimal value: %lX", hexadecimalval);		  
			}
			else if(bnm==4)
			{
				 char hexNum[100];
				long int count=0;
				printf("Enter a hexadecimal number To Convert it into Binary : ");
				scanf("%s",hexNum);
				printf("\nBinary Number is : ");
				while(hexNum[count])
				{
					switch(hexNum[count])
					{
						case '0' : printf("0000");
							break;
						case '1' : printf("0001");
							break;
						case '2' : printf("0010");
							break;
						case '3' : printf("0011");
							break;
						case '4' : printf("0100");
							break;
						case '5' : printf("0101");
							break;
						case '6' : printf("0110");
							break;
						case '7' : printf("0111");
							break;
						case '8' : printf("1000");
							break;
						case '9' : printf("1001");
							break;
						case 'A' : printf("1010");
							break;
						case 'B' : printf("1011");
							break;
						case 'C' : printf("1100");
							break;
						case 'D' : printf("1101");
							break;
						case 'E' : printf("1110");
							break;
						case 'F' : printf("1111");
							break;
						case 'a' : printf("1010");
							break;
						case 'b' : printf("1011");
							break;
						case 'c' : printf("1100");
							break;
						case 'd' : printf("1101");
							break;
						case 'e' : printf("1110");
							break;
						case 'f' : printf("1111");
							break;
						default : printf("\nInvalid Entry, Please Try Again  %c",hexNum[count]);
					}
					count++;
				}		  
			}
			else if(bnm==5)
			{
				long long binr;
    			printf("Enter a binary number: ");
    			scanf("%lld", &binr);
    			int oct = 0, dec = 0, i = 0;
    			while (binr != 0) {
        			dec += (binr % 10) * pow(2, i);
        			++i;
        			binr /= 10;
    			}
    			i = 1;
    			while (dec != 0) {
        			oct += (dec % 8) * i;
        			dec /= 8;
        			i *= 10;
    				}
    			printf("%lld in binary = %d in octal", binr, oct);
			}
			else if(bnm==6)
			{
				int oct;
    			printf("Enter an octal number: ");
    			scanf("%d", &oct);
    			int dec = 0, i = 0;
    			long long bin = 0;
    			while (oct != 0) {
        			dec += (oct % 10) * pow(8, i);
      				  ++i;
       				 oct /= 10;
    				}
   				 i = 1;
    			while (dec != 0) {
       			 bin += (dec % 2) * i;
        			dec /= 2;
       			 i *= 10;
   					 }
   				printf("%d in octal = %lld in binary", oct, bin);
			}
			else
			{
				printf("\nEnter a valid number.....\nchoose again.....");
				goto ibincon;
			}
	}
	else{
		printf("you entered a wrong value.....\nchoose again.....");
		goto ocon;
		}
		goto z;
			break;
		case 10:
			J:
			printf("\nSOME ALGEBRAIC FUNCTIONS:");
			saf:
			printf("\n1)SUM UPTO N NUMBERS\n2)SUM OF SQUARES UPTO N NUMBERS\n3)SUM OF CUBES UPTO N NUMBERS");
			printf("\n4)SUM OF POWER 4 TERMS UPTO N NUMBERS\n5)FACTORIAL OF A NUMBER\n6)NPR\n7)NCR\n8)GCD\n9)LCM");
			printf("\nEnter the operation to be performed from the menu given above:");
			scanf("%d",&opt);
			int af,lim,Sum=0,ar,an,max;
			if(opt>=1&&opt<=4)
			{
				printf("\nEnter the value of n upto which term summation takes place be:");
				scanf("%d",&lim);
				for(i=1;i<=lim;i++)
				{    
      				Sum=Sum+pow(i,opt);    
  				}
				printf("\nSum of n^%d terms upto %d terms:%d",opt,lim,Sum);
			}
			else if(opt==5)
			{
				printf("\nEnter the number whose factorial should be found:");
				scanf("%d",&af);
				printf("\nThe factorial of %d is:%d",af,factorial(af));
			}
			else if(opt==6||opt==7)
			{
				printf("\nEnter the value of n:");
				scanf("%d",&an);
				printf("\nEnter the value of r:");
				scanf("%d",&ar);
				if(opt==6)
				{
					printf("\nThe value of %dP%d is:%d",an,ar,factorial(an)/factorial(an-ar));
				}
				if(opt==7)
				{
					printf("\nThe value of %dC%d is:%d",an,ar,factorial(an)/(factorial(an-ar)*factorial(ar)));
				}
			}
			else if(opt==8||opt==9)
			{
			int n1,n2,i,gcd;
    		printf("Enter two positive integers(n1 and n2) resectively: ");
    		gcdp:
    		scanf("%d %d", &n1, &n2);
    		if(n1<0||n2<0)
    		{
    			printf("\nenter a positive number......");
    			goto gcdp;
			}
    		if(opt==8)
			{
    		for(i=1; i <= n1 && i <= n2; ++i)
    		{
        		if(n1%i==0 && n2%i==0)
            		gcd = i;
    		}
    		printf("G.C.D of %d and %d is %d", n1, n2, gcd);	
    		}
    		if(opt==9)
    		{
    			max = (n1 > n2) ? n1 : n2;
    			while (1) {
        			if (max % n1 == 0 && max % n2 == 0) {
            			printf("The LCM of %d and %d is %d.", n1, n2, max);
            			break;
        			}
        			++max;
    			}
			}
			}
			else
			{
				printf("\nEnter a valid number.....\nchoose again.....");
				goto saf;
			}
			goto z;
			break;
		case 11:
			K:
    		printf("\nBASIC GEOMETRY");
    		ogn:
    		printf("\n1)Plane figures\n2)3-d figures");
    		int gn,g1,g2;
    		float atr,gs,si1,si2,si3,let,gl,gb,gr,ga,gc,gd,gh,gbase;
    		printf("\nEnter an option from above :");
    		scanf("%d",&gn);
    		if(gn==1)
    		{
        		printf("PLANE FIGURES");
        		ig1:
        		printf("\n1)Square/Rhombus\n2)Rectangle/parallelogram\n3)Circle\n4)Triangle\n5)Traphesium");
        		printf("\nEnter an option from above :");
        		scanf("%d",&g1);
        		if(g1==1)
        		{
            		printf("\n....Square....");
            		printf("\nEnter the side of the square:");
            		scanf("%f",&gs);
            		printf("\nArea of square:%f",gs*gs);
            		printf("\nPerimeter of square:%f",4*gs);
        		}
        		else if(g1==2)
        		{
            		printf("\n.... rectangle...");
            		printf("\nEnter the length and breadth of the rectangle:");
            		scanf("%f%f",&gl,&gb);
            		printf("\nArea of rectangle:%f",gl*gb);
            		printf("\nPerimeter of rectangle :%f",2*(gl+gb));
        		}
        		else if(g1==3)
        		{
            		printf("\n.... circle.......");
            		printf("\nEnter the radius of the circle:");
            		scanf("%f",&gr);
            		printf("\nArea of circle :%f",pi*gr*gr);
            		printf("\n circumference of circle :%f",2*pi*gr);
        		}
        		else if(g1==4)
        		{
            		printf("\n.... Triangle......");
            		printf("\n1)by sides\n2)by base and height");
            		itro:
            		if(atr==1)
            		{
            			printf("\nEnter the sides of the triangle:");
            			scanf("%f%f%f",si1,si2,si3);
            			let=(si1+si2+si3)/3;
            			printf("\nArea of the triangle:%f",sqrt(let*(let-si1)*(let-si2)*(let-si3)));
            			printf("\nperimeter of the triangle:",let*3);
					}
            		else if(atr==2)
            		{
            		printf("\nEnter the base and height of a traingle :");
            		scanf("%f%f",&gbase,&gh);
            		printf("\nArea of triangle :%f",0.5*gbase*gh);
        			}
        			else
        			{
        				printf("\nEnter a valid number.....\nchoose again.....");
        				goto itro;
					}
        		}
        		else if(g1==5)
        		{
            		printf("\n....trapezium......");
            		printf("\nEnter the height and four sides(starting with base and upper side respectively) of the trapezium:");
            		scanf("%f%f%f%f%f",&gh,&ga,&gb,&gc,&gd);
            		printf("\nArea of trapezium :%f",0.5*gh*(ga+gb));
            		printf("\n perimeter of a trapezium  :%f",ga+gb+gc+gd);
        		}
        		else
        		{
        			printf("\nEnter a valid number....\nchoose again.....");
        			goto ig1;
				}
    		}
    		else if(gn==2)
    		{
		            printf("3-D FIGURES");
		            ig2:
		            printf("\n1)Cube\n2)Cuboid\n3)Cone\n4)Square pyramid\n5)Prism\n6)Cylinder\n7)Sphere\n8)Hemisphere");
		            printf("\nEnter an option from above :");
		            scanf("%d",&g2);
		            if(g2==1)
        		    {
                		printf("\n.......CUBE......");
        		        printf("\nEnter the length of the edge:");
                		scanf("%f",&gc);
                		printf("\nLateral surface area of the cube:%f",4*gc*gc);
            		    printf("\nTotal surface area of the Cube:%f",6*gc*gc);
    		            printf("\nVolume of the cube:%f",gc*gc*gc);
            		}
            		else if(g2==2)
            		{
                		printf("\n.......CUBOID.....");
                		printf("\nEnter the length, breadth and height:");
                		scanf("%f%f%f",&gl,&gb,&gh);
                		printf("\nLateral surface area of the Cuboid:%f",2*gh*(ga+gb));
                		printf("\nTotal surface area of the Cuboid:%f",2*(gl*gb+gb*gh+gh*gl));
                		printf("\nVolume of the cuboid:%f ",gl*gb*gh);
            		}
            		else if(g2==3)
            		{
                		printf("\n.......CONE.....");
                		printf("\nEnter the length,height and radius:");
                		scanf("%f%f%f",&gl,&gh,&gr);
                		printf("\nLateral surface area of the Cone:%f",pi*gr*gl);
                		printf("\nTotal surface area of the Cone:%f",(pi*gr*gr)+(pi*gr*gl));
                		printf("\nVolume of the cone: %f",1/3*(pi*gr*gr*gh));
		            }
		        	else if(g2==4)
		        	{
		            printf("\n.......square pyramid.....");
		            printf("\nEnter the length of the base edge and height:");
		            scanf("%f%f",&gl,&gh);
		            printf("\nLateral surface area of the square pyramid:%f",2*gl*sqrt((gl*gl/4)+gh*gh));
		            printf("\nTotal surface area of the square pyramid :%f",(gl*gl)+(2*gl)*sqrt((gl*gl/4)+gh*gh));
		            printf("\nVolume of the square cuboid: %f",1/3*(gl*gl+gh));
		        	}
		        	else if(g2==5)
		        	{
		            printf("\n......PRISM.....");
		            printf("\nEnter the length of base edge and height of the prism and length of the prism respectively:");
		            scanf("%f%f%f",&gb,&gh,&gl);
		            printf("\nEnter the lengths of the triangular edges:");
		            scanf("%f%f%f",&ga,&gb,&gc);
		            printf("\nLateral surface area of the prism :%f",(ga+gb+gc)*gl);
		            printf("\nTotal surface area of the prism :%f",(ga+gb+gc)*gl+0.5*gb*gh*gl);
		            printf("\nVolume of the prism: %f",0.5*gb*gh*gl);
        			}
        			else if(g2==6)
        			{
            		printf("\n.......CYLINDER....");
            		printf("\nEnter the radius and height:");
            		scanf("%f%f",&gr,&gh);
            		printf("\nLateral surface area of the cylinder :%f",2*pi*gr*gh);
            		printf("\nTotal surface area of the cylinder :%f",2*pi*gr*(gh+gr));
            		printf("\nVolume of the cylinder:%f",pi*gr*gr*gh);
        			}
        			else if(g2==7)
        			{
            		printf("\n........SPHERE.......");
            		printf("\nEnter the radius:");
            		scanf("%f",&gr);
            		printf("\nLateral surface area of the sphere :%f",4*pi*gr*gr);
            		printf("\nTotal surface area of the sphere :%f",4*pi*gr*gr);
            		printf("\nVolume of the sphere:%f ",(4/3)*(pi*gr*gr*gr));
        			}
        			else if(g2==8)
        			{
        			printf("\n........HEMI-SPHERE.......");
            		printf("\nEnter the radius:");
            		scanf("%f",&gr);
            		printf("\nLateral surface area of the hemi-sphere :%f",2*pi*gr*gr);
            		printf("\nTotal surface area of the hemi-sphere :%f",3*pi*gr*gr);
            		printf("\nVolume of the hemi-sphere:%f",(2/3)*(pi*gr*gr*gr));
        			}	
        			else
        			{
        				printf("\nEnter a valid number......\nchoose again.....");
        				goto ig2;
					}
    		}
    		else
    		{
        		printf("\nEnter a valid number....\nchoose again.....\n");
        		goto ogn;
    		}
    		goto z;
			break;
		case 12:
			printf("\nSCIENTIFIC CONSTANTS:");
			printf("\nproton mass 1.673289821*10^-27kg mp");
			printf("\nneutron mass 1.67492747121*10^-27kg mp");
			printf("\nelectron mass 9.1093835611*10^-31");
			printf("\nmuon mass 1.88353159448*10^-28 kg");
			printf("\nbohr radius 5.291772106712*10^-11 m");
			printf("\nplanck constant 6.62607004081*10^-34 js");
			printf("\nnuclear magneton 5.05078369931*10^-27JTA^-1");
			printf("\nbohr magneton 277400999457*10^-24 JTA-1");
			printf("\nplanck constant. rationalised 1.05457180013*10^-34 Js");
			printf("\nfine-structure constant 0.007297352566417 ");
			printf("\nclassical electron radius 2.817940322719*10^-15 m");
			printf("\ncompton wavelength 2.426310236711*10^-12 m");
			printf("\nproton gyromagnetic ratio 2.67522190018*10^8 SA-1T^-1");
			printf("\nproton compton wavelength 1.3214098539661*10^-15 m");
			printf("\nneutron compton wavelength 1.3195909048188*10^-15m");
			printf("\nrydberg constant 1.09737315680865*10^7 m^-1");
			printf("\natomic mass constant 1.6605390402*10^-27 kg");
			printf("\nElectron magnetic moment  -9.28476462057*10^-24 JT^-1");
			printf("\nproton magnetic moment  1.410606787397*10^-26JT^-1");
			printf("\nneutron magnetic moment  -9.662365023*10^-27 JT^-1");
			printf("\nmuon magnetic moment  -4.490448261*10^-26 JT^-1");
			printf("\nfaraday constant  96485.332895 Cmol^-1");
			printf("\nelementary charge  1.602176620898*10^-19 C");
			printf("\navogadro constant  6.02214085774*10^23 mol^-1");
			printf("\nboltsmann constant 1.3806485279*10^-23 JK^-1 f");
			printf("\nmolar volume of ideal 0.02271094713m^3mol^-1");
			printf("\nmolar gas constant  8.3144598484848jmol^-1k^-1");
			printf("\nspeed of light in vaccum 2.99792458*10^8 ms^-1");
			printf("\nfirst radiation constant    3.74177179046*10^-16   Wm^2");
			printf("\nsecond radiation constant   0.014387773683 mk");
			printf("\nstefan-boltsmann constant    5.67036713*10^-8 Wm^-2K^-4");
			printf("\nelectric constant       8.854187817*10^-12  Fm^-1");
			printf("\nmagnetic constant       1.25663706*10^-6 NA^-2");
			printf("\nmagnetic flux  quantum   2.06783383113*10^-15 Wb");
			printf("\nstandard accelaration of gravity    9.80665ms^-2");
			printf("\nconductance  quantum    7.748091731018*10^-5  Wb");
			printf("\ncharateristic impedence of vaccum    376.730313461 ohms");
			printf("\ncelsius temperature       273.15 K");
			printf("\nnewtonian constant of gravitation  6.6740831*10^-11 m^3kg^-1s^-2");
			printf("\nstandard atmosphere     101325.0 Pa");
			goto z;
			break;
		default:
			printf("\nENTER A VALID NUMBER.....\nchoose again\n");
			goto head;
	}
	z:
		printf("\n\n\nEnter any option from below list to continue:");
		printf("\npress'1'for END\npress'2' for HOME PAGE\npress'3'forBACK");
				scanf("%d",&des);
				if(des==1)
				{	
					printf("\n\t\t------****END****------");
					exit (0);
				}
				else if(des==2)
				{
					goto head;
				}
				else if(des==3)
				{
					if(n==1)
					{
						goto A;
					}
					else if(n==2)
					{
						goto B;
					}
					else if(n==3)
					{
						goto C;
					}
					else if(n==4)
					{
						goto D;
					}
					else if(n==5)
					{
						goto E;
					}
					else if(n==6)
					{
						goto F;
					}
					else if(n==7)
					{
						goto G;
					}
					else if(n==8)
					{
						goto H;
					}
					else if(n==9)
					{
						goto I;
					}
					else if(n==10)
					{
						goto J;
					}
					else if(n==11)
					{
						goto K;
					}
					else if(n==12)
					{
						goto head;
					}
				}	
				else
					printf("\n Enter a valid number.....");
}
