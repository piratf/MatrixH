##简易矩阵运算库

实现了部分矩阵运算函数，目标是实现尽可能快，参考了部分开源实现代码。
所有方法放在一个头文件中实现，只需要 include 头文件即可，便于使用。

编译指令：

```
make
```
或
```
g++ -std=c++11 -O3
```
<!-- 对于 `Matrixomp.h` 编译指令：g++ -std=c++11 -O3 --openmp 使 openmp 生效 -->

方法详见头文件类声明部分，详细注释在实现部分

#####需要注意的部分

> * warning: 取模运算符 % 被重载为 右除，便于使用。（不是矩阵取模
> * warning: 可以使用泛型实例化，但使用 `complex<T>` 类型实例化矩阵时由于 `toDouble` 函数的原因无法使用除法等相关函数。计划通过继承实现类的拆分。
> * /= 运算由于在 除法部分 涉及到数值的类型转换，暂时没有优化，其他赋值运算符有优化，比运算再赋值快。