module Spice


export furnsh,
       unload,
       pxform,
       utc2et,
       str2et,
       vsep,
       spkpos,
       spkpos!,
       timout,
       reclat


function getSpiceLib()
 sharedLib = ""
 iFile = open("../.userSettings.conf")
 while !eof(iFile)
   line = readline(iFile)
   if contains(line, "spicelib:")
     sharedLib = string(bytestring(split(line, "spicelib:")[2][1:end-1]))
     return sharedLib
   end
 end
end

function furnsh(KernelFile::AbstractString)
  ccall((:furnsh_c, sharedLib),Void,(Ptr{Cchar},),KernelFile)
end

function pxform(from, to, et)
  rotMat = zeros(Float64, 3,3)
  ccall((:pxform_c, sharedLib), Void, (Ptr{Cchar}, Ptr{Cchar}, Cdouble, Ptr{Cdouble}),
        from, to, et, rotMat)
  return rotMat

end

function reclat(a::Vector)
  r = Ref{Cdouble}(0.0)
  lon = Ref{Cdouble}(0.0)
  lat = Ref{Cdouble}(0.0)
  ccall((:reclat_c, sharedLib), Void, (Ptr{Cdouble}, Ptr{Cdouble},
                     Ptr{Cdouble}, Ptr{Cdouble}), a, r, lon, lat)
  return r[], lon[], lat[]
end


function str2et(s::AbstractString)
  et = Ref{Cdouble}(0.0)
  ccall((:utc2et_c, sharedLib), Float64, (Ptr{Cchar}, Ptr{Cdouble}), s, et)
  return et[]
end

function scs2e(sc::Int64, sclkch::ASCIIString)
  et = Ref{Cdouble}(0.0)
  ccall((:scs2e_c, sharedLib), Void, (Int64, Ptr{Cchar}, Ptr{Cdouble}),
        sc, sclkch, et)
  return et[]
end

function spkpos(target::ASCIIString, et::Float64, ref::ASCIIString,
                abcorr::ASCIIString, obs::ASCIIString)
  pos = zeros(Float64, 3)
  lt = Ref{Cdouble}(0.0)
  ccall((:spkpos_c, sharedLib), Void, (Ptr{Cchar}, Float64,
        Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cdouble}, Ptr{Cdouble}),
        target, et, ref, abcorr, obs, pos, lt)
  return pos, lt[]
end

function spkpos!(target::ASCIIString, et::Float64, ref::ASCIIString,
                abcorr::ASCIIString, obs::ASCIIString, pos)
  lt = Ref{Cdouble}(0.0)
  ccall((:spkpos_c, sharedLib), Void, (Ptr{Cchar}, Float64,
        Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cdouble}, Ptr{Cdouble}),
        target, et, ref, abcorr, obs, pos, lt)
  nothing
end

function timout(et::Float64, pictur::ASCIIString, lenOutput::Int64)
  s = "............................."
  ccall((:timout_c, sharedLib), Void, (Float64, Ptr{Cchar}, Int64, Ptr{Cchar}),
        et, pictur, lenOutput, s)
  return s
end

function unload(KernelFile::ASCIIString)
  ccall((:furnsh_c, sharedLib),Void,(Ptr{Cchar},),KernelFile)
end

function utc2et(s::ASCIIString)
  et = Ref{Cdouble}(0.0)
  ccall((:utc2et_c, sharedLib), Float64, (Ptr{Cchar}, Ptr{Cdouble}), s, et)
  return et[]
end

function vsep(a::Vector, b::Vector)
  ccall((:vsep_c, sharedLib), Float64, (Ptr{Float64}, Ptr{Float64}), a, b)
end

const sharedLib = getSpiceLib()

end
