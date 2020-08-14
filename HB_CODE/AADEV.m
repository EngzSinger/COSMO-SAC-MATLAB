function ad = AADEV(a)

ave = mean(a)*ones(length(a),1);
dev = abs(a-ave);
ad = sum(dev)/mean(a);
end
