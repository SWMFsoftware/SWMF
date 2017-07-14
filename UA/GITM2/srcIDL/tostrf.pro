function tostrf,value
  return, strcompress(string(float(value)),/remove_all)
end
