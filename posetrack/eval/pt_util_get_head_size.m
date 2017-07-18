function headSize = pt_util_get_head_size(rect)

SC_BIAS = 1;
headSize = SC_BIAS*norm([rect.x2 rect.y2] - [rect.x1 rect.y1]);

end