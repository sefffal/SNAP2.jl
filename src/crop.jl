
## Crop masters to size of lights (centred)
# Assumes A is larger than B, and we are cropping it down.
function crop_A_to_B(A, B)
    if size(A,1) <= size(B,1) && size(A,2) <= size(B,2)
        return A
    end
    xr_A = axes(A,1)
    yr_A = axes(A,2)
    xr_A = xr_A .- mean(xr_A)
    yr_A = yr_A .- mean(yr_A)

    xr_B = axes(B,1)
    yr_B = axes(B,2)
    xr_B = xr_B .- mean(xr_B)
    yr_B = yr_B .- mean(yr_B)

    m1 = xr_B[begin] .<= xr_A .<= xr_B[end]
    m2 = yr_B[begin] .<= yr_A .<= yr_B[end]
    m1 = findfirst(m1):findlast(m1)
    m2 = findfirst(m2):findlast(m2)
    A_crop = copyheader(A, A[m1,m2])
    return A_crop
end
export crop_A_to_B