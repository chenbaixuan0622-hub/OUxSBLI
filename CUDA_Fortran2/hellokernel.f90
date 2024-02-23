program main
    use kernel
    use cudafor
    implicit none

    integer :: a
    integer,device :: dev_a        ! デバイスメモリの割付（デバイス変数の宣言）

    a = 1
    dev_a = a                      ! GPUメモリへ値をコピー
    call doublify<<<1,1>>>(dev_a)  ! カーネル呼び出し
    a = dev_a                      ! GPUメモリから結果をコピー
    print *,a
end program main