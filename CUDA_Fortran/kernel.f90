module kernel
    implicit none
    contains
    attributes(global) subroutine doublify(a) ! attributes(global)を付ける
        integer,intent(inout) :: a
        a = a*2
    end subroutine doublify
end module kernel