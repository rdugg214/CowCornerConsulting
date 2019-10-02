function helper_plotcmap(X,Z,MAT);
figure
surf(X,Z,MAT)
view(2)
caxis([-10 -5])
colorbar

end