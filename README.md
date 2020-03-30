# GGB2MMA

相比 Mathematica, GeoGebra 在教师中受众更广，（免费，开源）

本文提供 GGB2MMA 的代码实现，由于 GGB 的函数重载现象很普遍，所以每个不同的功能都用单独写一行来区分

如果需要在函数之间找对应，只需 `Ctrl`+`F` 在网页搜索

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
|   9   |  atan(x) 或 arctan(x)  |                Arc tangent，反正切函数（返回值在 -π/2 与π/2 之间）                 |                   ArcTan[x]                   |
|  10   | atand(x) 或 arctand(x) |                     Arc tangent，反正切 (tan-1)（以度为单位）                      |               ArcTan[x Degree]                |
|  11   | atanh(x) 或 arctanh(x) |                       Antihyperbolic tangent，反双曲正切函数                       |                  ArcTanh[x]                   |
|  12   |       atan2(y,x)       |                    Arctangent，反正切函数（返回值在 -π与π之间）                    | ArcTan[x,y]  (\*这里 x,y 顺序和前面是反的 \*) |
|  13   |       beta(a,b)        |                             Β(a,b)，贝塔函数（β函数）                              |                   Beta[a,b]                   |
|  14   |      beta(a,b,x)       |                              Β(x;a,b)，不完全贝塔函数                              |                  Beta[x,a,b]                  |
|  15   | betaRegularized(a,b,x) |                           I(x;a,b)，正则化不完全贝塔函数                           |           BetaRegularized[x, a, b]            |
|  16   |        cbrt(x)         |                            Cubicroot，三次方根、立方根                             |                  CubeRoot[x]                  |
|  17   |        ceil(x)         |                 Least integer greater than or equal，“向上取整”，                  |                  Ceiling[x]                   |
|  18   |      conjugate(x)      |                                Conjugate，共轭函数                                 |                 Conjugate[x]                  |
|  19   |         cos(x)         |                                      余弦函数                                      |                    Cos[x]                     |
|  20   |   cosec(x) 或 csc(x)   |                   Cosecant，余割函数  cosec∠A=c/a（斜边 / 对边）                   |                    Csc[x]                     |
|  21   |  cosech(x) 或 csch(x)  |                         Hyperbolic cosecant，双曲余割函数                          |                    Csch[x]                    |
|  22   |        cosh(x)         |                           Hyperbolic cosin，双曲余弦函数                           |                    Cosh[x]                    |
|  23   |     cosIntegral(x)     |                             Cosine Integral，余弦积分                              |                CosIntegral[x]                 |
|  24   |         cot(x)         |                   Cotangent，余切函数   cot∠A=c/a（邻边 / 对边）                   |                    Cot[x]                     |
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
|  45   |         sec(x)         |                     Secant，正割函数  sec∠A=c/b（斜边 / 邻边）                     |                    Sec[x]                     |
|  46   |        sech(x)         |                          Hyperbolic secant，双曲正割函数                           |                    Sech[x]                    |
|  47   |   sgn(x) 或 sign(x)    |                Sign，符号（x 为正数返回 1，负数返回 -1，零返回 0）                 |                    Sign[x]                    |
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

|   ID   |          GeoGebra          |                     名称                     |                         Mathematica                          |
| :----: | :------------------------: | :------------------------------------------: | :----------------------------------------------------------: |
| 2.1.1  |      Side(\<Quadric>)      | 侧面 (《二次曲面》) , 创建有限二次曲面的侧面 | Graphics3D[{CapForm[None], Tube[{{0, 0, 0}, {0, 0, 3}}, 1]}] |
| 2.1.2  |     PerpendicularPlane     |                   垂直平面                   |                   InfinitePlane[p,{v1,v2}]                   |
| 2.1.8  | Sphere(\<Point>,\<Radius>) |                     球面                     |                  Sphere[<Point>,\<Radius>]                   |
| 2.1.9  |          Surface           |                     曲面                     |                       ParametricPlot3D                       |
| 2.1.11 |           Volume           |                     体积                     |                            Volume                            |
| 2.1.16 |          Cylinder          |                     圆柱                     |                           Cylinder                           |
| 2.1.1  |            Cone            |                     圆锥                     |                             Cone                             |
| 2.1.16 |          Cylinder          |                     圆柱                     |                           Cylinder                           |

</center>

## 2.2 指令

<center>

|  ID   | GeoGebra  |       名称       |                                          Mathematica                                           |
| :---: | :-------: | :--------------: | :--------------------------------------------------------------------------------------------: |
|       |  ToPoint  |     转换为点     |                                              ReIm                                              |
|       | ToComplex |    转换为复数    |                                          {x,y}.{1,I}I                                          |
|       |  ToPolar  | 转换为极坐标形式 | 如果复数，CoordinateTransformData["Cartesian" -> "Polar", "Mapping", ReIm@x], 向量就去掉 ReIm@ |

</center>

## 2.5 代数

|  ID   |                    GeoGebra                    |    名称    |                             Mathematica                              |
| :---: | :--------------------------------------------: | :--------: | :------------------------------------------------------------------: |
|       |                     Cross                      |    叉积    |                                Cross                                 |
|       |                    Division                    |    除法    |                          QuotientRemainder                           |
|       |                    Division                    | 多项式除法 |                     PolynomialQuotientRemainder                      |
|       |                      Dot                       |            |                                 Dot                                  |
|       |               CommonDenominator                |            |                  PolynomialLCM@Denominator /@ list                   |
|       |                   NextPrime                    |            |                              NextPrime                               |
|       |                    Simplify                    |            |                  Simplify      或     FullSimplify                   |
|       |                   Solutions                    |            |                            Values@*Solve                             |
|       |                     NSolve                     |            |                                NSolve                                |
|       | NSolve(\<Equation>,\<Variable=starting value>) |            |                               FindRoot                               |
|       |                   NSolutions                   |            |                            Values@*NSolve                            |
|       |                 CompleteSquare                 |            | /. a_. x_^2 + b_. x_ + c_ :> a (x + b/(2 a))^2 + (4 a c - b^2)/(4 a) |
|       |                 PreviousPrime                  |            |                           NextPrime[x,-1]                            |
|       |                     Solve                      |            |                                Solve                                 |
|       |                      Mod                       |            |                                 Mod                                  |
|       |                      Mod                       |            |                            PolynomialMod                             |
|       |                      Div                       |            |                               Quotient                               |
|       |                      Div                       |            |                          PolynomialQuotient                          |
|       |                    IFactor                     |            | Factor[x^2 - x - 1, Extension -> All]     (\*只能在复数域上分解 \*)  |
|       |                    IsPrime                     |            |                                PrimeQ                                |
|       |                     Factor                     |            |                                Factor                                |
|       |                    Divisors                    |            |                          DivisorSigma[0,x]                           |
|       |                  DivisorsSum                   |            |             DivisorSum[x, # &]   或  DivisorSigma[1, x]              |
|       |                  DivisorsList                  |            |                               Divisors                               |
|       |                   RightSide                    |            |                             formula[[2]]                             |
|       |                   RightSide                    |            |                  list/. Equal -> List //#[[;;,-1]]&                  |
|       |                     Expand                     |            |                                Expand                                |
|       |                  PrimeFactors                  |            |            ConstantArray @@@ FactorInteger[x] // Flatten             |
|       |                     ToBase                     |            |                         IntegerDigits[23,8]                          |
|       |                    FromBase                    |            |                              FromDigits                              |
|       |                      GCD                       |            |                                 GCD                                  |
|       |                      GCD                       |            |                            PolynomialGCD                             |
|       |                      Max                       |            |                                 Max                                  |
|       |                      Max                       |            |     list // Transpose // ConstantArray @@@ # & // Flatten // Max     |
|       |                      Max                       |            |                             FindMaximum                              |
|       |                      LCM                       |            |                                 LCM                                  |
|       |                      LCM                       |            |                            PolynomialLCM                             |
|       |                      Min                       |            |                                 Min                                  |
|       |                      Min                       |            |     list // Transpose // ConstantArray @@@ # & // Flatten // Min     |
|       |                      Min                       |            |                             FindMaximum                              |
|       |                    LeftSide                    |            |                             formula[[1]]                             |
|       |                    LeftSide                    |            |                  list/. Equal -> List //#[[;;,1]]&                   |

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



## 2.7 Function.函数与微积分


<center>

|                   GeoGebra                    |            Mathematica            |
| :-------------------------------------------: | :-------------------------------: |
|               PartialFractions                |               Apart               |
|             ParametricDerivative              |                 D                 |
|                    Degree                     |             Exponent              |
|                  Derivative                   |                                   |
|                Iteration.迭代                 |               Nest                |
|            IterationList.迭代列表             |             NestList              |
|               NIntegral.定积分                |            NIntegrate             |
|             Polynomial.多项式函数             |              Expand               |
|         Polynomial(\<List of Points>)         |      InterpolatingPolynomial      |
| Polynomial(\<Function>)；多项式函数(<函数>)。 |      //Expand//Collect[#,x]&      |
|     NInvert(\<Function>)；反函数(<函数>)      |          InverseFunction          |
|                  Denominator                  |            Denominator            |
|                   Numerator                   |             Numerator             |
|                  ComplexRoot                  |               Solve               |
|                InflectionPoint                |       D[#,x]&//Solve[#==0]&       |
|                   Normalize                   |             Normalize             |
|                   Function                    |        InterpolationOrder         |
|                   Function                    |           Plot  Plot3D            |
|                   Integral                    |             Integrate             |
|                IntegralBetween                |   Integrate[f[x]-g[x],{x,a,b}]    |
|                     Limit                     |               Limit               |
|                   Extremum                    |             Maximize              |
|                   Asymptote                   |               放弃                |
|                   SolveODE                    |              DSolve               |
|                   SolveODE                    |              NDSolve              |
|                   NSolveODE                   |              NDSolve              |
|                 RectangleSum                  |
|                     Root                      |               Solve               |
|                     Root                      |             FindRoot              |
|                   RootList                    |           {#,0}&/@list            |
|                 PathParameter                 |   点值 （Mathematica 无此概念）   |
|               OsculatingCircle                |              密切圆               |
|                SVD.奇异值分解                 |    SingularValueDecomposition     |
|                Curvature.曲率                 |           ArcCurvature            |
|                CurvatureVector                |
|         Curve（CurveCartesian）.曲线          |  ParametricPlot ParametricPlot3D  |
|            TrigCombine.三角式合并             |            TrigReduce             |
|                 TrigSimplify                  |           FullSimplify            |
|                  TrigExpand                   |            TrigExpand             |
|                   UpperSum                    |
|                  NDerivative                  |                 D                 |
|               TaylorPolynomial                |              Series               |
|            TrapezoidalSum.梯形法则            |
|                 Coefficients                  |          CoefficientList          |
|                 LowerSum.下和                 |
|               SlopeField.斜率场               |            VectorPlot             |
|                Spline.样条曲线                |           BSplineCurve            |
|                    Factors                    |              Factor               |
|                    Factors                    |           FactorInteger           |
|            ImplicitCurve.隐式曲线             |
|          ImplicitDerivative.隐式微分          |
|               LimitAbove.右极限               | Limit[, Direction -> "FromAbove"] |
|                 LeftSum.左和                  |                                   |
|               LimitBelow.左极限               | Limit[, Direction -> "FromBelow"] |
|                                               |

</center>

2.8 Geometry.几何

|                                           GeoGebra                                            |                     Mathematica                     |
| :-------------------------------------------------------------------------------------------: | :-------------------------------------------------: |
|                                          Radius.半径                                          |              ArcLength[Circle[]]/(2Pi)              |
|                                         envelope.包络                                         |                       无对应                        |
|                                        Difference.差异                                        |                  RegionDifference                   |
|                                          Length.长度                                          |                      ArcLength                      |
|                                          Length.长度                                          |                       Length                        |
|                           PerpendicularLine（OrthogonalLine）.垂线                            |                       无对应                        |
|                             Vertex(\<Conic>)；顶点(<圆锥曲线>)。                              |                       无对应                        |
|                            Vertex(\<Inequality>)；顶点(<不等式>)。                            |                       无对应                        |
|                             Vertex(\<Polygon>)；顶点(<多边形>)。                              |                       无对应                        |
|                         Polygon(\<List of Points>)；多边形(<点列>)。                          |                       Polygon                       |
|                      Direction(\<Line>)；方向向量(<直线 , 射线,线段>)。                       |                       无对应                        |
|                       RigidPolygon(\<Polygon>)；刚体多边形(<多边形>)。                        |                       无对应                        |
|       Locus(\<Point Creating Locus Line Q>,\<Point P>)；轨迹(<构造轨迹的点>,<控制点>)。       |                       无对应                        |
|                       Locus(\<f(x,y)>,\<Point>)；轨迹(\<f(x,y)>,<点>)。                       |                       无对应                        |
|                          LocusEquation(\<Locus>)；轨迹方程(<轨迹>)。                          |                       无对应                        |
|     LocusEquation(\<Boolean Expression>,\<Free Point>)；轨迹方程(<布尔表达式>,<自由点>)。     |                       Boolean                       |
|                                           Arc.弧线                                            |                       无对应                        |
|                                        CrossRatio.交比                                        |                       无对应                        |
|                                        Intersect.交点                                         |                        Solve                        |
|                                    AngleBisector.角平分线                                     |                       无对应                        |
|                                         Distance.距离                                         |                   RegionDistance                    |
|                Angle(\<Object>)；角度(<对象 圆锥曲线,向量,点,数字,多边形>)）。                |                     AngleVector                     |
|                                           Area.面积                                           |                        Area                         |
|                            Point(\<List>)；描点(<有序数组列表>)。                             |                        Point                        |
|                                            PointIn                                            |                     RandomPoint                     |
|                         InteriorAngles(\<Polygon>)；内角(<多边形>)。                          |                       无对应                        |
|                      Tangent(\<Point>,\<Conic>)；切线(<点>,<圆锥曲线>)。                      |                       无对应                        |
|               ClosestPointRegion(\<Region>,\<Point>)；区域内最近点(<区域>,<点>)               |                    RegionNearest                    |
|          Cubic(\<Point>,\<Point>,\<Point>,\<Number>)；三次曲线(<点 1>,<点 2>,<点 3>           |                       无对应                        |
|                                 CircumcircularSector.三点扇形                                 |                       无对应                        |
|                                  CircumcircularArc.三点圆弧                                   |                       无对应                        |
|                                    TriangleCurve.三角曲线                                     |                       无对应                        |
| TriangleCenter(\<Point>,\<Point>,\<Point>,\<Number>)；三角形中心(<点 1>,<点 2>,<点 3>,<数字>) |                                                     |
|                                     Trilinear.三线坐标点                                      |                       无对应                        |
|             Sector(\<Conic>,\<Point>,\<Point>)；扇形(<圆或椭圆>,<点 1>,<点 2>)。              |                       无对应                        |
|                       Ray(\<Start Point>,\<Point>)；射线(<起点>,<点>)。                       |                       无对应                        |
|               ArePerpendicular(\<Line>,\<Line>)；是否垂直(<直线 1>,<直线 2>)。                |                       无对应                        |
|        AreConcurrent(\<Line>,\<Line>,\<Line>)；是否共点(<直线 1>,<直线 2>,<直线 3>)。         |                       无对应                        |
|          AreCollinear(\<Point>,\<Point>,\<Point>)，是否共线(<点 1>,<点 2>,<点 3>)。           |                       无对应                        |
|                                     AreConcyclic.是否共圆                                     |
|                  AreParallel(\<Line>,\<Line>)；是否平行(<直线 1>,<直线 2>。                   |
|            AreCongruent(\<Object>,\<Object>)；是否全等(<几何对象1>,<几何对象 2>)。            |
|             AreEqual(\<Object>,\<Object>)；是否相等(<几何对象 1>,<几何对象 2>)。              |                        Equal                        |
|                   IsTangent(\<Line>,\<Conic>)；是否相切(<直线>,<圆锥曲线>)                    |
|                       Segment(\<Point>,\<Point>)；线段(<点 1>,<点 2>)。                       |                        Line                         |
|                                    IntersectPath.相交路径                                     |
|                          Slope(\<Line>)；斜率(<直线, 射线, 线段>)。                           |
|                            Centroid(\<Polygon>)；形心(<多边形>)。                             |                   RegionCentroid                    |
|         CircularArc(\<Midpoint>,\<Point A>,\<Point B>)；圆弧(<圆心>,<点 1>,<点 2>)。          |
|       CircularSector(\<Midpoint>,\<Point A>,\<Point B>)；圆扇形(<圆心>,<点 1>,<点 2>)。       |
|                          Circumference(Conic)；圆周长(<圆锥曲线>)。                           |                      ArcLength                      |
|                          Polyline(\<List of Points>)；折线(<点列>)。                          |                        Line                         |
|                      Prove(\<Boolean Expression>)；证明(<布尔表达式>)。                       |
|                 ProveDetails(\<Boolean Expression>)；证明过程(<布尔表达式>)。                 |
|                        Line(\<Point>,\<Point>)；直线(<点 1>,<点 2>)。                         |                    InfiniteLine                     |
|                          PerpendicularBisector(LineBisector).中垂线                           |                PerpendicularBisector                |
|                                         Midpoint.中点                                         |                       无对应                        |
|           Barycenter(\<List of Points>,\<List of Weights>)；重心(<点列>,<权重列表>)           | //MapThread[ConstantArray, #]&//Flatten[#,1]&//Mean |
|                            Perimeter(\<Polygon>)；周长(<多边形>)。                            |                      Perimeter                      |
|                     ClosestPoint(\<Path>,\<Point>)；最近点(<路径>,<点>)。                     |                    RegionNearest                    |
|                 Reflect(\<Object>,\<Point>)；对称(<几何对象>,<对称中心点>)。                  |                 ReflectionTransform                 |
|                   Translate(\<Object>,\<Vector>)；平移(<几何对象>,<向量>)。                   |                      Translate                      |
|          Shear(\<Object>,\<Line>,\<Ratio>)；切变(<几何对象>,<直线,射线,线段>,<比>)。          |                  ShearingTransform                  |
|                    Stretch(\<Object>,\<Vector>)；伸缩(<几何对象>,<向量>)。                    |
|                                     Dilate(Enlarge).位似                                      |
|                                            Rotate                                             |                       Rotate                        |

## 2.10 ScriptingCommands.脚本指令

|                                       GeoGebra                                       |  Mathematica  |
| :----------------------------------------------------------------------------------: | :-----------: |
|                          Button(\<Caption>)；按钮("<标题>")                          |    Button     |
|                        PlaySound(\<URL>)；播放声音(<网址>)。                         |     Play      |
|                        PlaySound(\<URL>)；播放声音(<网址>)。                         |     Sound     |
|                                 ExportImage.导出图片                                 |    Export     |
|                     ZoomIn(\<Scale Factor>)；放大(<缩放因子>)。                      | Ctrl+鼠标滚轮 |
|                                    SetValue.赋值                                     |      Set      |
|                                   Checkbox.复选框                                    |   Checkbox    |
|                  CopyFreeObject(\<Object>)；复制自由对象(<对象>)。                   |
|     AttachCopyToView(\<Object>,\<View 0,1,2>)；附加副本(<对象>,<视图 0 ,1,2>)。      |
|                         UpdateConstruction( )；更新作图( )。                         |
|                                  Turtle()；海龟()。                                  |
|               TurtleBack(\<Turtle>,\<Distance>)；后退(<海龟>,<路程>)。               |
|                                    Slider.滑动条                                     |    Slider     |
|        ParseToFunction(\<Function>,\<String>)；解析为函数(<函数>,<字符串>)。         | ToExpression  |
|                                ParseToNumber.解析为数                                | ToExpression  |
|                                 StartRecord.开始记录                                 |
|                        TurtleDown(\<Turtle>)；落笔(<海龟>)。                         |
|                        Pan(\<x>,\<y>)；平移视图(\<x>,\<y>)。                         |
|                               StartAnimation.启动动画                                |
|                                StartLogging.启动日志                                 |
|             TurtleForward(\<Turtle>,\<Distance>)；前进(<海龟>,<路程>)。              |
|                          Delete(\<Object>)；删除(<对象>)。                           |
|                           SetBackgroundColor.设置背景颜色                            |
|                                SetDecoration.设置标记                                |
|                              SetLabelMode.设置标签模式                               |
|            SetCaption(\<Object>,\<Text>)；设置标题(<对象>,"<标题文本>")。            |
|                                SetPointSize.设置点径                                 |   PointSize   |
|                                SetPointStyle.设置点型                                |  PointLegend  |
|                             SetDynamicColor.设置动态颜色                             |
|     SetFixed(\<Object>,\<true,false>)；设置设置对象锁定(<对象>,\<true,false>)。      |
|                               SetPerspective.设置格局                                |
|                                  SetTrace.设置跟踪                                   |
|                           SetTooltipMode.设置工具提示模式                            |    Tooltip    |
|                              SetActiveView.设置活动视图                              |
|                             SetVisibleInView.设置可见性                              |
|                            SetViewDirection.设置视图方向                             |   ViewPoint   |
|                                 SetFilling.设置填充                                  |
|                                  SetLayer.设置图层                                   |
| SetLevelOfDetail(\<Surface>,\<Level of Detail>)；设置细节级别(<曲面>,<细节级别 0,1>) |
|                        SetConditionToShowObject.设置显示条件                         |
|                              SetLineThickness.设置线径                               |   Thickness   |
|                                SetLineStyle.设置线型                                 |    Dashed     |
|                                  SetColor.设置颜色                                   |
|                       SetSeed(\<Integer>)；设置种子(<整数>)。                        |  SeedRandom   |
|                                SetSpinSpeed.设置转速                                 |
|                                  SetCoords.设置坐标                                  |
|                             SetAxesRatio.设置坐标轴比例                              |  AspectRatio  |
|                             InputBox（Textfield）.输入框                             |  InputField   |
|                                DataFunction.数据函数                                 |
|                                     ZoomOut.缩小                                     |
|                                    TurtleUp.抬笔                                     |
|                                 StopLogging.停止日志                                 |
|                                   GetTime.系统时间                                   | TimeObject[]  |
|                                  ShowLabel.显示标签                                  |
|                       ShowLayer(\<Number>)；显示图层(<数值>)。                       |
|                              ShowGrid( )；显示网格( )。                              |   GridLines   |
|                                 ShowAxes.显示坐标轴                                  |     Axes      |
|                                  SelectObjects.选择                                  |
|                                  HideLayer.隐藏图层                                  |
|                                   TurtleRight.右转                                   |
|                        ReadText(\<Text>)；阅读文本("<文本>")                         |
|                             RunClickScript.运行单击脚本                              |
|                             RunUpdateScript.运行更新脚本                             |
|                    Execute(\<List of Texts>)；执行(<文本列表>)。                     | ToExpression  |
|    CenterView(\<Center Point>)；中心定位(<视图中心设置坐标(x,y) , 视图中心点>)。     |  AxesOrigin   |
|                                     Repeat.重复                                      |
|                                    Rename.重命名                                     |
|                                   TurtleLeft.左转                                    |

## 2.11 DiscreteMath.离散数学


|                           GeoGebra                           |   Mathematica    |
| :----------------------------------------------------------: | :--------------: |
|            DelaunayTriangulation.Delaunay 三角网             |   DelaunayMesh   |
|                      Voronoi.Voronoi 图                      |   VoronoiMesh    |
|                 TravelingSalesman.旅行商问题                 | FindShortestTour |
|        ConvexHull(\<List of Points>)；凸包(<点列>)。         |  ConvexHullMesh  |
|                  ShortestDistance.最短距离                   | FindShortestPath |
| MinimumSpanningTree(\<List of Points>)；最小生成树(<点列>)。 | FindSpanningTree |



## 2.12 列表

|   GeoGebra    |                      Mathematica                       |
| :-----------: | :----------------------------------------------------: |
|    Flatten    |                        Flatten                         |
|     Union     |                   Union RegionUnion                    |
|    Insert     |                         Insert                         |
|    Product    |                      Times@@list                       |
|    Product    |             Block[{Plus=Times},Total@list]             |
|    Product    |                        Product                         |
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
|     Take      |                          Take                          |
|     Take      |                       StringTake                       |
|   Sequence    |                         Table                          |
|   Sequence    |                         Range                          |
|      Zip      |     Map[f,list]         (\*区别是 mma 有纯函数 \*)     |
|    Element    |                          [[]]                          |
|    Append     |                         Append                         |
|     Last      |                          Last                          |
|     First     |                         First                          |

## 2.13 Logical.逻辑


<center>

|  GeoGebra  | Mathematica  |
| :--------: | :----------: |
|     If     |      If      |
| IsInteger  |   IntegerQ   |
| IsDefined  |  Definition  |
| IsInRegion | RegionMember |
|  CountIf   |    Count     |
|   KeepIf   |    Select    |

## 2.14 Statistics.统计

|    GeoGebra    |                                       Mathematica                                        |
| :------------: | :--------------------------------------------------------------------------------------: |
|      Mad       |                                      MeanDeviation                                       |
|      Mean      |                                           Mean                                           |
|    SigmaXX     | list//#[[;;,1]]&//#^2&//Total    或  list//#[[;;,1]]&// #.#& 或  list//First/@#& // #.#& |
|    SigmaXY     |                              list// Times @@@ # & // Total                               |
|    SigmaXY     |                                        list//#.#&                                        |
|    SigmaYY     |                                  list//Last/@#& // #.#&                                  |
|    Spearma     |                                       SpearmanRho                                        |
|  Stdev.Sstdev  |                                   list//Variance//Sqrt                                   |
| Stdevp.Stdevp  |                                    StandardDeviation                                     |
|      Sxx       |                          list// (#.# - (Total@# )^2)/Length@# &                          |
|      Sxx       |                    list// First/@#&//(#.# - (Total@# )^2/Length@#) &                     |
|      Sxy       |     list// Total@(Times @@@ #) - (Total@(First /@ #))*(Total@(Last /@ #))/Length@# &     |
|      Sxy       |             list// #1.#2 - (Total@#1*Total@#2)/Length@# &[#[[1]], #[[2]]] &              |
|      Syy       |                   list// Last /@ # & // (#.# - (Total@#)^2/Length@#) &                   |
|     TTest      |                                          TTest                                           |
|   Percentile   |                                         Quantile                                         |
|       SD       |                                    StandardDeviation                                     |
|   ZMeanTest    |                            ZTest     (\*ZTest 系列不清楚 \*)                             |
|       Q1       |                           Quartiles[{1, 3, 4, 2, 5, 6}]//First                           |
|       Q3       |                           Quartiles[{1, 3, 4, 2, 5, 6}]//Last                            |
|  HarmonicMean  |                                       HarmonicMean                                       |
|     FitLog     |       NonlinearModelFit[{{E, 1}, {E^2, 4}}, {a Log[x] + b}, {a, b},   x] // Normal       |
|    FitPoly     |            Fit[{{-1, -1}, {0, 1}, {1, 1}, {2, 5}}, x^Range[0, 3], x] // Chop             |
|    Variance    |                                         Variance                                         |
|     ANOVA      |                                          不清楚                                          |
|      SDX       |  横坐标标准差     {{1,1},{2,2},{3,1},{3,3},{4,2},{3,-1}} //First/@#&//StandardDeviation  |
|     MeanX      |                                  list//First/@#&//Mean                                   |
| GeometricMean  |                                      GeometricMean                                       |
| RootMeanSquare |                                      RootMeanSquare                                      |
| ChiSquaredTest |                                   PearsonChiSquareTest                                   |
|    RSquare     |                                          不清楚                                          |
|  FitLogistic   |               NonlinearModelFit[list, {a/(1 + b E^(-c x))}, {a, b, c}, x]                |
|     FitPow     |                       NonlinearModelFit[list, {a x^b}, {a, b}, x]                        |
|      Fit       |                                           Fit                                            |
|    FitLineX    |                                    Fit[list,{1,x},x]                                     |
|    FitLineY    |                                          不清楚                                          |
|  TTestPaired   |                                       PairedTTest                                        |
|   FitGrowth    |                     NonlinearModelFit[list, {a b ^( x)}, {a, b}, x]                      |
|    Shuffle     |                                       RandomSample                                       |
|   Covariance   |                                        Covariance                                        |
|     Sample     |                                       RandomSample                                       |
|  FitImplicit   |                                   拟合隐函数,没想出来                                    |
|     FitSin     |                NonlinearModelFit[list, {a+ b Sin[c x+d]}, {a, b,c,d}, x]                 |
|     FitExp     |                       NonlinearModelFit[list, {E ^(a x)}, {a}, x]                        |
|     Median     |                                          Median                                          |
|      Mode      |                                        Commonest                                         |
|      Sum       |                                          Total                                           |
|      Sum       |                                        StringJoin                                        |
|      SDY       |  纵坐标标准差     {{1,1},{2,2},{3,1},{3,3},{4,2},{3,-1}} //Last/@#&//StandardDeviation   |
|     MeanY      |                                   list//Last/@#&//Mean                                   |

</center>

## 2.15 Chart.图表


|                                GeoGebra                                 |        Mathematica        |
| :---------------------------------------------------------------------: | :-----------------------: |
|              StickGraph(\<List of Points>)；棒图(<点列>)。              | ListPlot[Filling -> Axis] |
|  ResidualPlot(\<List of Points>,\<Function>)；残差图(<点列>,<函数>)。   | list/.{a_,b_}:>{a,b-f@a}  |
|         DotPlot(\<List of Raw Data>)；点阵图(<原始数据列表>)。          |                           |
|             StepGraph(\<List of Points>)；阶梯图(<点列>)。              |
|                   StemPlot(\<List>)；茎叶图(<列表>)。                   |
|                         ContingencyTable.列联表                         |
|                          FrequencyTable.频数表                          |
|                       FrequencyPolygon.频数多边形                       |
|                             BarChart.条形图                             |         BarChart          |
|                             BoxPlot.箱线图                              |
| NormalQuantilePlot(\<List of Raw Data>)；正态分位数图(<原始数据列表>)。 |
|                            Histogram.直方图                             |         Histogram         |
|                         HistogramRight.直方图右                         |


## 2.16 Text.文本

|                             GeoGebra                              |     Mathematica      |
| :---------------------------------------------------------------: | :------------------: |
| TableText(\<List>,\<List>,...)；表格文本(<列表 1>,<列表 2>,...)。 |      TableForm       |
|            FractionText(\<Number>)；分数文本(<数字>)。            |     Rationalize      |
|              SurdText(\<Number>)；根式文本(<数值>)。              |   RootApproximant    |
|            FormulaText(\<Object>)；公式文本(<对象>)。             |        MaTeX         |
|          ScientificText(\<Number>)；科学计数法(<数字>)。          |    ScientificForm    |
|          ContinuedFraction(\<Number>)；连分式(<数字>)。           |  ContinuedFraction   |
|            VerticalText(\<Text>)；竖排文本("<文本>")。            | Text[]~Rotate~(Pi/2) |
|       UnicodeToText(\<List of Integers>)；统一码转换为文本        |  FromCharacterCode   |
|                 UnicodeToLetter.统一码转换为字母                  |  FromCharacterCode   |
|                  Text(\<Object>)；文本(<对象>)。                  |         Text         |
|      TextToUnicode("\<Text>")；文本转换为统一码("<文本>")。       |   ToCharacterCode    |
|               Ordinal(\<Integer>)；序数(<自然数>)。               |       ToString       |
|  RotateText(\<Text>,\<Angle>)；旋转文本("<文本>",<角度,弧度>)。   |     Text//Rotate     |
|    LetterToUnicode("\<Letter>")；字母转换为统一码("<字母>")。     |   ToCharacterCode    |


## 2.17 Vector&Matrix.向量与矩阵

|                          GeoGebra                           |     Mathematica      |
| :---------------------------------------------------------: | :------------------: |
| UnitPerpendicularVector(\<Line>)；单位法向量(<直线,射线>)。 |      Normalize       |
|           Identity(\<Number>)；单位矩阵(<数值>)。           |    IdentityMatrix    |
|          UnitVector(\<Vector>)；单位向量(<向量>)。          |      Normalize       |
|       PerpendicularVector(\<Line>)；法向量(<直线>)。        | {a,b}/.{a,b}:>{-b,a} |
|  ReducedRowEchelonForm(\<Matrix>)；简化行梯阵式(<矩阵>)。   |      RowReduce       |
|          MatrixRank(\<Matrix>)；矩阵的秩(<矩阵>)。          |      MatrixRank      |
|              Invert(\<Matrix>)；逆反(<矩阵>)。              |       Inverse        |
|             Invert(\<Function>)；逆反(<函数>)。             |   InverseFunction    |
|        Dimension(\<Object>)；维度(<点\|向量\|矩阵>)         |      Dimensions      |
|        Vector(\<Point>)；向量(<终点(原点为起点)>)。         |       VectorQ        |
|          Determinant(\<Matrix>)；行列式(<矩阵>)。           |         Det          |
| ApplyMatrix(\<Matrix>,\<Object>)；应用矩阵(<矩阵>,<对象>)。 |         Dot          |
|            Transpose(\<Matrix>)；转置(<矩阵>)。             |      Transpose       |

## 2.18 Optimization.优化指令

|                                  GeoGebra                                   | Mathematica |
| :-------------------------------------------------------------------------: | :---------: |
| Maximize(\<Dependent number>,\<Free number>)；最大值点(<因变数>,<滑动条>)。 |  NMaximize  |
|  Minimize(\<Dependent number>,\<Free number>)；最小值点(<因变量>,<参数>)。  |  NMinimize  |

## 2.19 Conic.圆锥曲线

|                               GeoGebra                               | Mathematica |
| :------------------------------------------------------------------: | :---------: |
|          LinearEccentricity(\<Conic>)；半焦距(<圆锥曲线>)。          |
|         Semicircle(\<Point>,\<Point>)；半圆(<点 1>,<点 2>)。         |
|        SemiMinorAxisLength(\<Conic>)；副半轴长(<圆锥曲线>)。         |
|               MinorAxis(\<Conic>)；副轴(<圆锥曲线>)。                |
| ConjugateDiameter(\<Vector>,\<Conic>)；共轭直径(<向量>,<圆锥曲线>)。 |
|          Polar(\<Point>,\<Conic>)；极线(<点>,<圆锥曲线>)。           |
|               Parameter(\<Parabola>)；焦参数(<抛物线>)               |
|                 Focus(\<Conic>)；焦点(<圆锥曲线>)。                  |
|             Eccentricity(\<Conic>)；离心率(<圆锥曲线>)。             |
| Incircle(\<Point>,\<Point>,\<Point>)；内切圆(<点 1>,<点 2>,<点 3>)。 |  Insphere   |
|         Parabola(\<Point>,\<Line>)；抛物线(<焦点>,<准线>)。          |
|                           Hyperbola.双曲线                           |
|                             Ellipse.椭圆                             |
|                             Circle.圆周                              |   Circle    |
|                            Conic.圆锥曲线                            |
|                         Center(Centre).中心                          |
|                              Axes.轴线                               |
|           SemiMajorAxisLength（FirstAxisLength）.主半轴长            |
|                     MajorAxis（FirstAxis）.主轴                      |
|                            Directrix.准线                            |

## 2.20 CASSpecific.运算区专属指令

|                                      GeoGebra                                      |                       Mathematica                        |
| :--------------------------------------------------------------------------------: | :------------------------------------------------------: |
|                      MixedNumber(<Number>)；带分数(<数值>)。                       |
| GroebnerDegRevLex(<List of Polynomials>)；分次反字典序 Groebner 基(<多项式列表>)。 |                      GroebnerBasis                       |
|   GroebnerLexDeg(<List of Polynomials>)；分次字典序 Groebner 基(<多项式列表>)。    |                      GroebnerBasis                       |
|                        CSolve(<Equation>)；复数解(<方程>)。                        |                          Solve                           |
|                     CSolutions(<Equation>)；复数解集(<方程>)。                     |                      Solve//Values                       |
|                 CFactor(<Expression>)；复数域因式分解(<表达式>)。                  |                  Factor[Extension->All]                  |
|               CIFactor(<Expression>)；复无理数域因式分解(<表达式>)。               |                  Factor[Extension->All]                  |
|            SolveCubic(<Cubic Polynomial>)；解三次多项式(<三次多项式>)。            |                          Solve                           |
|                     Numeric(<Expression>)；近似数(<表达式>)。                      |                            N                             |
|                    Laplace(<Function>)；拉普拉斯变换(<函数>)。                     |                     LaplaceTransform                     |
|                InverseLaplace(<Function>)；拉普拉斯逆变换(<函数>)。                |                 InverseLaplaceTransform                  |
|   Substitute(<Expression>,<from>,<to>)；替换(<表达式>,<被替换对象>,<替换对象>)。   |                            /.                            |
|  Eliminate(<List of Polynomials>,<List of Variables>)；消元(<多项式集>,<变量集>)   |                        Eliminate                         |
|                      Rationalize(<Number>)；有理化(<数值>)。                       |                       Rationalize                        |
|             ToExponential(<Complex Number>)；转换为指数形式(<复数>)。              | With[{n = Abs[1 + I], a = Arg[1 + I]}, Defer[n E^(I a)]] |
|       GroebnerLex(<List of Polynomials>)；字典序 Groebner 基(<多项式列表>)。       |                      GroebnerBasis                       |