function tobinary(filename)
	mf = matfile(filename);
	props = properties(mf);
	suj = props{2};
	
	pre = mf.(suj)(1,1);
	pre = pre{1}';
	pre = pre(:,384:-1:1);
	fid=fopen("../0_data/pre.bin","w"); fwrite(fid,pre,'double'); fclose(fid);
	clear pre;
	
	mov_pre = mf.(suj)(1,2);
	mov_pre = mov_pre{1};
	dlmwrite("../0_data/mov_pre.dat",mov_pre,' ');
	clear mov_pre;
	
	post = mf.(suj)(2,1);
	post = post{1}';
	post = post(:,384:-1:1);
	fid=fopen("../0_data/post.bin","w"); fwrite(fid,post,'double'); fclose(fid);
	clear post;
	
	mov_post = mf.(suj)(2,2);
	mov_post = mov_post{1};
	dlmwrite("../0_data/mov_post.dat",mov_post,' ');
	clear mov_post;
end

