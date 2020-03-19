
相比 Mathematica, GeoGebra 在教师中受众更广，（免费，开源）

## 1.2 函数
<center>

|  ID   |        GeoGebra        |                                        名称                                        |                  Mathematica                  |
| :---: | :--------------------: | :--------------------------------------------------------------------------------: | :-------------------------------------------: |
|   1   |         abs(x)         |                               Absolute value，绝对值                               |                    Abs[x]                     |
|   2   |  acos(x) 或 arccos(x)  |                               Arc cosine，反余弦函数                               |                   ArcCos[x]                   |
|   3   | acosd(x) 或 arccosd(x) |                      Arc cosine，反余弦 (cos-1)（以度为单位）                      |               ArcCos[x Degree]                |
|   4   | acosh(x) 或 arccosh(x) |                       Antihyperbolic cosine，反双曲余弦函数                        |                  ArcCosh[x]                   |
|   5   |         arg(x)         |                                Argument，复数的幅角                                |                    Arg[z]                     |
|   6   |  asin(x) 或 arcsin(x)  |                                Arc sine，反正弦函数                                |                   ArcSin[x]                   |
|   7   | asind(x) 或 arcsind(x) |                       Arc sine，反正弦 (sin-1)（以度为单位）                       |               ArcSin[x Degree]                |
|   8   | asinh(x) 或 arcsinh(x) |                        Antihyperbolic sine，反双曲正弦函数                         |                  ArcSinh[x]                   |
|   9   |  atan(x) 或 arctan(x)  |                 Arc tangent，反正切函数（返回值在-π/2 与π/2 之间）                 |                   ArcTan[x]                   |
|  10   | atand(x) 或 arctand(x) |                     Arc tangent，反正切 (tan-1)（以度为单位）                      |               ArcTan[x Degree]                |
|  11   | atanh(x) 或 arctanh(x) |                       Antihyperbolic tangent，反双曲正切函数                       |                  ArcTanh[x]                   |
|  12   |       atan2(y,x)       |                    Arctangent，反正切函数（返回值在-π与π之间）                     | ArcTan[x,y]  (\*这里 x,y 顺序和前面是反的 \*) |
|  13   |       beta(a,b)        |                             Β(a,b)，贝塔函数（β函数）                              |                   Beta[a,b]                   |
|  14   |      beta(a,b,x)       |                              Β(x;a,b)，不完全贝塔函数                              |                  Beta[x,a,b]                  |
|  15   | betaRegularized(a,b,x) |                           I(x;a,b)，正则化不完全贝塔函数                           |           BetaRegularized[x, a, b]            |
|  16   |        cbrt(x)         |                            Cubicroot，三次方根、立方根                             |                  CubeRoot[x]                  |
|  17   |        ceil(x)         |                 Least integer greater than or equal，“向上取整”，                  |                  Ceiling[x]                   |
|  18   |      conjugate(x)      |                                Conjugate，共轭函数                                 |                 Conjugate[x]                  |
|  19   |         cos(x)         |                                      余弦函数                                      |                    Cos[x]                     |
|  20   |   cosec(x) 或 csc(x)   |                    Cosecant，余割函数  cosec∠A=c/a（斜边/对边）                    |                    Csc[x]                     |
|  21   |  cosech(x) 或 csch(x)  |                         Hyperbolic cosecant，双曲余割函数                          |                    Csch[x]                    |
|  22   |        cosh(x)         |                           Hyperbolic cosin，双曲余弦函数                           |                    Cosh[x]                    |
|  23   |     cosIntegral(x)     |                             Cosine Integral，余弦积分                              |                CosIntegral[x]                 |
|  24   |         cot(x)         |                    Cotangent，余切函数   cot∠A=c/a（邻边/对边）                    |                    Cot[x]                     |
|  25   |        coth(x)         |                         Hyperbolic cotangent，双曲余切函数                         |                    Coth[x]                    |
|  26   |    exp(x) 或 $e^x$     |                           Exponential function，指数函数                           |                    Exp[x]                     |
|  27   |         erf(x)         |                       Gaussian Error Function，高斯误差函数                        |                    Erf[x]                     |
|  28   |     expIntegral(x)     |                           Exponential Integral，指数积分                           |                ExpIntegral[x]                 |
|  29   |   fractionalPart(x)    |                             Fractional Part，小数函数                              |               FractionalPart[x]               |
|  30   |        floor(x)        |                  Greatest integer less than or equal，“向下取整”                   |                   Floor[x]                    |
|  31   |        gamma(x)        |                            Γ(x)，伽玛函数（Gamma 函数）                            |                   Gamma[x]                    |
|  32   |       gamma(a,x)       |                           γ(a,x)，（低阶）不完全伽玛函数                           |                  Gamma[a,x]                   |
|  33   | gammaRegularized(a,x)  |                  P(a,x)=γ(a,x)/Γ(a)，（低阶）正则化不完全伽玛函数                  |             GammaRegularized[a,x]             |
|  34   |      imaginary(x)      |                        Imaginary，虚值函数（参见实值函数）                         |                     Im[x]                     |
|  35   |         ld(x)          |                      Logarithm to base2，以 2 为底的对数函数                       |              Log2[x] 或 Log[2,x]              |
|  36   |         lg(x)          |                     Logarithm to base10，以 10 为底的对数函数                      |             Log10[x] 或 Log[10,x]             |
|  37   |    ln(x) 或 log(x)     |                            Logarithm，自然对数，底数为ℯ                            |                    Log[x]                     |
|  38   |        log(b,x)        |                    Logarithm of x to baseb，以 b 为底 x 的对数                     |                   Log[b,x]                    |
|  39   |       nroot(x,n)       |                                  求 x 的 n 次方根                                  |            Surd[x, n]  或 x^(1/n)             |
|  40   |     polygamma(m,x)     |                          Polygamma function，多项伽玛函数                          |                 PolyGamma[x]                  |
|  41   |         psi(x)         |            Digamma function，Digamma 函数（伽玛函数的对数的导数），Ψ(x)            |                PolyGamma[0, x]                |
|  42   |        real(x)         |                              实值函数，复数的实部数值                              |                     Re[x]                     |
|  43   |      random(a,b)       |                 Random number between a and b，a 与 b 之间的随机数                 |               RandomReal[a, b]                |
|  44   |        round(x)        |                                  Round，四舍五入                                   |                   Round[x]                    |
|  45   |         sec(x)         |                      Secant，正割函数  sec∠A=c/b（斜边/邻边）                      |                    Sec[x]                     |
|  46   |        sech(x)         |                          Hyperbolic secant，双曲正割函数                           |                    Sech[x]                    |
|  47   |   sgn(x) 或 sign(x)    |                 Sign，符号（x 为正数返回 1，负数返回-1，零返回 0）                 |                    Sign[x]                    |
|  48   |         sin(x)         |                                   Sine，正弦函数                                   |                    Sin[x]                     |
|  49   |        sinh(x)         |                           Hyperbolic sine，双曲正弦函数                            |                    Sinh[x]                    |
|  50   |     sinIntegral(x)     | Sine Integral，正弦积分       缩写 Si(z)，Si(z)=∫sin(t)/tdt, 区间 (0,z) 上的定积分 |                SinIntegral[x]                 |
|  51   |        sqrt(x)         |                           Square root，平方（二次方）根                            |                    Sqrt[x]                    |
|  52   |         tan(x)         |                                 Tangent，正切函数                                  |                    Tan[x]                     |
|  53   |        tanh(x)         |                          Hyperbolic tangent，双曲正切函数                          |                    Tanh[x]                    |
|  54   |        zeta(x)         |                                ζ(x)，黎曼 zeta 函数                                |                    Zeta[x]                    |

以下几个非严格的函数，是 GeoGebra 内部约定的计算或变量
|  ID   |                GeoGebra                 |                                 名称                                 | Mathematica |
| :---: | :-------------------------------------: | :------------------------------------------------------------------: | :---------: |
|  55   |                  x(x)                   | x-coordinate，点对象 x 的横坐标值           x(A) 返回点 A 的横坐标值 |   First@A   |
|  56   |                  xAxis                  |  x 轴                 非严格意义的函数，不需要变量，相当于直线 y=0   |    y==0     |
|  57   |                  y(x)                   | y-coordinate，点对象 x 的纵坐标值             Y(A) 返回点 A 的纵坐标 |   A[[2]]    |
|  58   |                  yAxis                  |  y 轴                 非严格意义的函数，不需要变量，相当于直线 x=0   |    x==0     |
|  59   | z(x):z-coordinate，点对象 x 的 z 坐标值 |                      Z(A) 返回点 A 的 z 坐标值                       |   Last@A    |
|  60   |                  zAxis                  |  z 轴                 非严格意义的函数，不需要变量，相当于 z 轴直线  | x==0&&y==0  |

</center>

另外，设有列表“表 1={1,2,3,a,D}”,“表 1(3)”返回列表的第三号元素“3”
表 1[[3]]

## 1.3 布尔运算

略

## 2.1 3D. 三维

> 几何作图部分没有直接对应，所以我就自己发挥了
<center>

|   ID   |         GeoGebra         |                    名称                    |                         Mathematica                          |
| :----: | :----------------------: | :----------------------------------------: | :----------------------------------------------------------: |
| 2.1.1  |     Side(<Quadric>)      | 侧面 (<二次曲面>) , 创建有限二次曲面的侧面 | Graphics3D[{CapForm[None], Tube[{{0, 0, 0}, {0, 0, 3}}, 1]}] |
| 2.1.2  |    PerpendicularPlane    |                  垂直平面                  |                   InfinitePlane[p,{v1,v2}]                   |
| 2.1.8  | Sphere(<Point>,<Radius>) |                    球面                    |                   Sphere[<Point>,<Radius>]                   |
| 2.1.9  |         Surface          |                    曲面                    |                       ParametricPlot3D                       |
| 2.1.11 |          Volume          |                    体积                    |                            Volume                            |
| 2.1.16 |         Cylinder         |                    圆柱                    |                           Cylinder                           |
| 2.1.1  |           Cone           |                    圆锥                    |                             Cone                             |
| 2.1.16 |         Cylinder         |                    圆柱                    |                           Cylinder                           |

</center>

## 2.2 指令

<center>

|  ID   | GeoGebra  |          名称           |                                          Mathematica                                           |
| :---: | :-------: | :---------------------: | :--------------------------------------------------------------------------------------------: |
|       |  ToPoint  |    转换为点 (<复数>)    |                                              ReIm                                              |
|       | ToComplex | 转换为复数 (<向量或点>) |                                          {x,y}.{1,I}I                                          |
|       |  ToPolar  |    转换为极坐标形式     | 如果复数，CoordinateTransformData["Cartesian" -> "Polar", "Mapping", ReIm@x], 向量就去掉 ReIm@ |

</center>

## 2.5 代数

|  ID   |                   GeoGebra                   | 名称  |                                   Mathematica                                    |
| :---: | :------------------------------------------: | :---: | :------------------------------------------------------------------------------: |
|       |                    Cross                     | 叉积  |                                      Cross                                       |
|       | Division(<Dividend Number>,<Divisor Number>) | 除法  |                 QuotientRemainder 或 PolynomialQuotientRemainder                 |
|       |                     Dot                      |       |                                       Dot                                        |
|       |              CommonDenominator               |       |                        PolynomialLCM@Denominator /@ list                         |
|       |                  NextPrime                   |       |                                    NextPrime                                     |
|       |                   Simplify                   |       |                        Simplify      或     FullSimplify                         |
|       |                  Solutions                   |       |                                  Values@*Solve                                   |
|       |                    NSolve                    |       |                                      NSolve                                      |
|       | NSolve(<Equation>,<Variable=starting value>) |       |                                     FindRoot                                     |
|       |                  NSolutions                  |       |                                  Values@*NSolve                                  |
|       |                CompleteSquare                |       |       /. a_. x_^2 + b_. x_ + c_ :> a (x + b/(2 a))^2 + (4 a c - b^2)/(4 a)       |
|       |                PreviousPrime                 |       |                                 NextPrime[x,-1]                                  |
|       |                    Solve                     |       |                                      Solve                                       |
|       |                     Mod                      |       |                               Mod 或 PolynomialMod                               |
|       |                     Div                      |       |                   Quotient         或      PolynomialQuotient                    |
|       |                   IFactor                    |       |         Factor[x^2 - x - 1, Extension -> All]     (*只能在复数域上分解*)         |
|       |                   IsPrime                    |       |                                      PrimeQ                                      |
|       |                    Factor                    |       |                                      Factor                                      |
|       |                   Divisors                   |       |                                DivisorSigma[0,x]                                 |
|       |                 DivisorsSum                  |       |                   DivisorSum[x, # &]   或  DivisorSigma[1, x]                    |
|       |                 DivisorsList                 |       |                                     Divisors                                     |
|       |                  RightSide                   |       |                formula[[2]] 或 list/. Equal -> List //#[[;;,-1]]&                |
|       |                    Expand                    |       |                                      Expand                                      |
|       |                 PrimeFactors                 |       |         FactorInteger  或 ConstantArray @@@ FactorInteger[x] // Flatten          |
|       |                    ToBase                    |       |                               IntegerDigits[23,8]                                |
|       |                   FromBase                   |       |                                    FromDigits                                    |
|       |                     GCD                      |       |                              GCD     PolynomialGCD                               |
|       |                     Max                      |       | Max   list // Transpose // ConstantArray @@@ # & // Flatten // Max  FindMaximum  |
|       |                     LCM                      |       |                               LCM   PolynomialLCM                                |
|       |                     Min                      |       | Min   list // Transpose // ConstantArray @@@ # & // Flatten // Min   FindMaximum |
|       |                   LeftSide                   |       |               formula[[1]]    或 list/. Equal -> List //#[[;;,1]]&               |
|       |                                              |       |                                                                                  |

## 2.6 概率

> PDF 给出概率密度函数，CDF 给出累积分布函数

|       GeoGebra       |                               Mathematica                                |
| :------------------: | :----------------------------------------------------------------------: |
|    FDistribution     |    FRatioDistribution[n, 10], x]   CDF[FRatioDistribution[n, 10], x]     |
|    TDistribution     | PDF[StudentTDistribution[\[Nu]], x]  CDF[StudentTDistribution[\[Nu]], x] |
|        Erlang        |   PDF[ErlangDistribution[k, .3], x]  CDF[ErlangDistribution[k, .3], x]   |
|      Bernoulli       |                          BernoulliDistribution                           |
|       Poisson        |                           PoissonDistribution                            |
|    RandomPoisson     |                    RandomVariate@*PoissonDistribution                    |
|    HyperGeometric    |                        HypergeometricDistribution                        |
|      LogNormal       |                          LogNormalDistribution                           |
|     BinomialDist     |                           BinomialDistribution                           |
|       Binomial       |                                 Binomial                                 |
|       Uniform        |                           UniformDistribution                            |
|    RandomUniform     |                    RandomVariate@*UniformDistribution                    |
|      ChiSquared      |                          ChiSquareDistribution                           |
|        Cauchy        |                            CauchyDistribution                            |
|    RandomDiscrete    |              RandomChoice[{0.2, 0.2, 0.6} -> {a, b, c}, 20]              |
|       Logistic       |                           LogisticDistribution                           |
| InverseFDistribution |                                                                          |

## 2.12 列表

|   GeoGebra    |                      Mathematica                       |
| :-----------: | :----------------------------------------------------: |
|    Flatten    |                        Flatten                         |
|     Union     |                   Union RegionUnion                    |
|    Insert     |                         Insert                         |
|    Product    | Times@@list   Block[{Plus=Times},Total@list]  Product  |
|    Product    | // Transpose // # /. {a_, b_} :> a^b & // Times @@ # & |
|    Product    |       // MapThread[#1^#2 &, #] & // Times @@ # &       |
|     Join      |                          Join                          |
|    Unique     |                         Union                          |
| Intersection  |                      Intersection                      |
|    Reverse    |                        Reverse                         |
|   Frequency   |   // Tally // (SortBy[#, First] &) // #[[;; , -1]] &   |
|    Remove     |                       Complement                       |
|     Sort      |                          Sort                          |
| RandomElement |                      RandomChoice                      |
|    IndexOf    |                     FirstPosition                      |
|     Take      |                  Take     StringTake                   |
|   Sequence    |                     Table   Range                      |
|      Zip      |     Map[f,list]         (\*区别是 mma 有纯函数 \*)     |
|    Element    |                          [[]]                          |
|    Append     |                         Append                         |
|     Last      |                          Last                          |
|     First     |                         First                          |
