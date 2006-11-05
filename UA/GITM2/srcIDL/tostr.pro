function tostr,value
  return, strcompress(string(long(value)),/remove_all)
end
