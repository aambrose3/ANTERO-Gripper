function g = initFcn(s,L)
    g = [-0.2*cos((pi/(L))*s);
         (pi/(5*L))*sin((pi/(L))*s)];
end