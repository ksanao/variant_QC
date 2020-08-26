# Setup history
export HISTSIZE=1000
export HISTIGNORE="&:[bf]g:exit:pwd:ls:history:vi"
export HISTCONTROL=ignoredups

# enable programmable completion
if [ -f /etc/bash_completion ]
then
    ./etc/bash_completion
fi

# customise bash prompt
export PS1='\e[0;33m\]\u\e[36m\]@\e[36m\]\h:\e[1;32m\]\w\e[0;32m\][$(git branch 2>/dev/null | grep "^*" | colrm 1 2)]$ '

# customise path
export "PATH=/home/src/pipeline:$PATH"

