# Notion Clone

A web-based collaborative workspace similar to Notion, built with React, Node.js, Express, and MongoDB.

## Features

- User authentication (register, login, logout)
- Create, edit, and delete documents
- Rich text editing with formatting options
- Real-time collaboration
- Document sharing
- Hierarchical document organization

## Prerequisites

- Node.js (v14 or higher)
- MongoDB
- npm or yarn

## Installation

1. Clone the repository:
```bash
git clone <repository-url>
cd notion-clone
```

2. Install backend dependencies:
```bash
npm install
```

3. Install frontend dependencies:
```bash
cd client
npm install
cd ..
```

4. Create a `.env` file in the root directory with the following variables:
```
MONGODB_URI=mongodb://localhost:27017/notion-clone
JWT_SECRET=your-secret-key
PORT=5000
```

## Running the Application

1. Start the MongoDB server

2. Start the backend server:
```bash
npm run dev
```

3. In a new terminal, start the frontend development server:
```bash
cd client
npm start
```

The application will be available at:
- Frontend: http://localhost:3000
- Backend API: http://localhost:5000

## Technologies Used

- Frontend:
  - React
  - Material-UI
  - Slate.js (rich text editor)
  - Socket.io-client
  - Axios

- Backend:
  - Node.js
  - Express
  - MongoDB
  - Mongoose
  - Socket.io
  - JWT for authentication

## Project Structure

```
notion-clone/
├── client/                 # React frontend
│   ├── public/
│   └── src/
│       ├── components/     # React components
│       ├── contexts/       # React contexts
│       └── App.js          # Main App component
├── models/                 # MongoDB models
├── routes/                 # API routes
├── middleware/             # Express middleware
├── server.js              # Backend entry point
└── package.json           # Project dependencies
```

## Contributing

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add some amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the LICENSE file for details.