function [ta] = total_area(V,F)

ta = 0.5*sum(doublearea(V,F));