module gauss_quad_mod

  use kinds_mod

  implicit none

  private

  public gaussian_legendre

contains

  subroutine gaussian_legendre(n, x, w)

    integer     , intent(in ) :: n
    real(r8), intent(out) :: x(n), w(n)

    real(16) fn(n), ak(n)
    real(16) m
    integer i, j

    j = 0 ! 赋值控制循环变量的初值           
    m = -1.000001 ! 设置计算域[-1，1] 的下限，即代替-1 
    do i = 1, 200000 ! 这个循环次数应该是由步长0.00001决 定,计算方法：200000=2/0.000001     
      if (legendre_p(m, n) * legendre_p(m+0.00001, n) < 0) then ! 从下限处开始往上逐步累加，
        ! 由步长0.00001说明最多求解10^5个解
        j = j + 1 ! 记录这是第几个解
        fn(j) = bisect(m, m + 0.00001, n)
        ! 调用二分法求解程序在分好的一小段上求解，
        ! 将解存储在fn(j)
        ak(j) = 2.0 / (n * dlegendre_p1(fn(j), n) * dlegendre_pn(fn(j), n)) ! 高斯点的权重
        x(j) = fn(j)
        w(j) = ak(j)
      end if
      m = m + 0.00001 ! 执行完一次判断m向前推进一步
    end do

  end subroutine gaussian_legendre

  real(16) function legendre_p(x, n) ! 定义Legendre函数
    
    real(16), intent(in) :: x
    integer, intent(in) :: n
      
    real(16) a(n) ! 代表n阶勒让德多项式
    integer i

    a(1) = x ! 1阶勒让德多项式
    a(2) = 1.5 * x**2 - 0.5 ! 2阶勒让德多项式
    do i = 3, n
      a(i) = (2 * i - 1) * x * a(i-1) / i - (i - 1) * a(i-2) / i
      ! 利用递推关系产生n阶勒让德多项式
    end Do
    legendre_p = a(n) ! 生成的n阶勒让德多项式   
  
  end function legendre_p
    
  real(16) function dlegendre_p1(x, n) ! 生成的n-1阶勒让德多项式
    
    real(16), intent(in) :: x
    integer, intent(in) :: n

    real(16) a(n) ! a(n-1)代表n-1阶勒让德多项式
    integer i
    
    a(1) = x
    a(2) = 1.5 * x**2 - 0.5
    do i = 3, n - 1
      a(i) = (2 * i - 1) * x * a(i-1) / i - (i-1) * a(i-2) / i
    end do
    dlegendre_p1 = a(n-1) 
  
  end function dlegendre_p1
    
  real(16) function dlegendre_pn(x, n) ! 生成n阶勒让德多项式的导数表达式
    
    real(16),intent(in) :: x
    integer,intent(in) :: n
      
    real(16) a(n)
    integer i

    a(1) = x
    a(2) = 1.5 * x**2 - 0.5
    do i = 3, n
      a(i) = (2 * i - 1) * x * a(i-1) / i - (i - 1) * a(i-2) / i
    end do
    dlegendre_pn = n * a(n-1) / (1 - x**2) - n * x * a(n) / (1 - x**2)

  end function dlegendre_pn

  real(16) function bisect(a_in, b_in, n) ! 二分法求解函数的解

    real(16), intent(in) :: a_in, b_in
    integer, intent(in) :: n

    real(16) a, b, c

    ! a, b是传递进来的划分好的有一个解存在的区间
    a = a_in
    b = b_in
    do
      c = (a + b) / 2.0
      if (legendre_p(c, n) * legendre_p(a, n) < 0) then
        b = c
      else
        a = c
      end if
      if (b - a < 1.e-16) exit 
    end do
    bisect = c

  end function bisect

end module gauss_quad_mod
