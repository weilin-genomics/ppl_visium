# 使用Python官方镜像作为基础镜像
FROM python:3.9-slim

# 设置工作目录
WORKDIR /app

# 安装系统依赖
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    gcc \
    python3-dev \
    && rm -rf /var/lib/apt/lists/*

# 复制requirements.txt
COPY requirements.txt .

# 升级pip
RUN pip install --no-cache-dir --upgrade pip

# 安装Python依赖
RUN pip install --no-cache-dir -r requirements.txt --verbose

# 复制应用程序文件
COPY Spatial_Visium.py .
COPY pages/ ./pages/
COPY Examples/ ./Examples/

# 设置Streamlit的文件上传大小限制为1GB
ENV STREAMLIT_SERVER_MAX_UPLOAD_SIZE=1000

# 暴露Streamlit默认端口
EXPOSE 1066

# 启动Streamlit应用
CMD ["streamlit", "run", "Spatial_Visium.py", "--server.port", "1066", "--server.address=0.0.0.0"] 