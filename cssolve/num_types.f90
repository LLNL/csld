!!<summary>Kind descriptors for the precision of integer and real values.</summary>
MODULE num_types
implicit none
private
public dp, sp, si, li
!!<member name="dp">Double precision kind for decimal precision of at least 15
!!digits, a decimal exponent range of at least 307.</member>
!!<member name="sp">Single precision kind for decimal precision of at least 6
!!digits, a decimal exponent range of at least 37.</member>
!!<member name="si">Short kind for integers that can represent all
!!values ranging from -10^1 (exclusive) to 10^1 (exclusive).</member>
!!<member name="li">Long kind for integers that can represent all
!!values ranging from -10^18 (exclusive) to 10^18 (exclusive).</member>
integer, parameter:: dp=selected_real_kind(15,307)
integer, parameter:: sp=selected_real_kind(6,37)
integer, parameter:: si=selected_int_kind(1)
integer, parameter:: li=selected_int_kind(18)
END MODULE
